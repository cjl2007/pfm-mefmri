#!/usr/bin/env python3
"""Template-anchored relabeling for existing PFM Infomap communities.

This is a comparison/visualization post-processing step. It does not alter the
Infomap community segmentation. For each density column, it reserves low label
IDs for ordered network templates and assigns each reserved ID to at most one
large community by spatial-prior overlap. Remaining communities keep distinct
high IDs.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import List, Sequence, Tuple

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.cifti2.cifti2_axes import ScalarAxis


def mat_string(x) -> str:
    arr = np.asarray(x)
    if arr.size == 1:
        x = arr.item()
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="ignore").strip()
    if isinstance(x, np.ndarray):
        vals = x.ravel()
        if vals.dtype.kind in ("U", "S"):
            return "".join(str(v) for v in vals).strip()
        if vals.size == 1:
            return mat_string(vals[0])
    s = str(x).strip()
    return s if s else "Network"


def load_priors(path: Path) -> Tuple[np.ndarray, List[str], np.ndarray]:
    mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    if "Priors" not in mat:
        raise ValueError(f"{path} does not contain a Priors struct")
    pri = mat["Priors"]
    if not hasattr(pri, "Spatial"):
        raise ValueError("Priors.mat must contain Priors.Spatial")
    spatial = np.asarray(pri.Spatial, dtype=np.float64)
    if spatial.ndim != 2:
        raise ValueError("Priors.Spatial must be 2D")
    n_net = spatial.shape[1]

    if hasattr(pri, "NetworkLabels"):
        raw = np.asarray(pri.NetworkLabels, dtype=object).ravel()
        labels = [mat_string(x) for x in raw[:n_net]]
    else:
        labels = []
    labels += [f"Network {i:02d}" for i in range(len(labels) + 1, n_net + 1)]

    if hasattr(pri, "NetworkColors"):
        colors = np.asarray(pri.NetworkColors, dtype=np.float64)
        colors = np.atleast_2d(colors)[:n_net, :3]
        if colors.size and float(np.nanmax(colors)) > 1.0:
            colors = colors / 255.0
    else:
        colors = np.zeros((0, 3), dtype=np.float64)
    if colors.shape[0] < n_net:
        rng = np.random.default_rng(17)
        pad = rng.uniform(0.25, 0.85, size=(n_net - colors.shape[0], 3))
        colors = np.vstack([colors, pad])
    return np.clip(np.nan_to_num(spatial, nan=0.0), 0.0, None), labels[:n_net], np.clip(colors[:n_net], 0.0, 1.0)


def cortical_grayordinates(axis) -> np.ndarray:
    idx: List[int] = []
    for name, slc, _ in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            stop = slc.stop if slc.stop is not None else axis.size
            idx.extend(range(slc.start, stop))
    return np.asarray(idx, dtype=np.int64)


def save_dscalar(ref_img: nib.Cifti2Image, values: np.ndarray, out_path: Path) -> None:
    axis = ref_img.header.get_axis(1)
    hdr = nib.Cifti2Header.from_axes((ScalarAxis(["labels"]), axis))
    data = np.asarray(values, dtype=np.float32)[None, :]
    nib.save(nib.Cifti2Image(data, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def label_color(value: int, base_colors: np.ndarray, n_net: int) -> Tuple[int, int, int]:
    if 1 <= value <= n_net:
        rgb = np.rint(base_colors[value - 1] * 255.0).astype(int)
        return int(rgb[0]), int(rgb[1]), int(rgb[2])
    rng = np.random.default_rng(int(value) + 1009)
    rgb = rng.integers(50, 225, size=3)
    return int(rgb[0]), int(rgb[1]), int(rgb[2])


def write_label_list(path: Path, rows: Sequence[dict], labels: Sequence[str], colors: np.ndarray) -> None:
    n_net = len(labels)
    seen = sorted({int(r["new_label_id"]) for r in rows if int(r["new_label_id"]) > 0})
    with path.open("w", encoding="utf-8") as f:
        for value in seen:
            if 1 <= value <= n_net:
                name = labels[value - 1]
            else:
                old = next((int(r["old_community_id"]) for r in rows if int(r["new_label_id"]) == value), value)
                name = f"Community_{old:04d}"
            r, g, b = label_color(value, colors, n_net)
            f.write(f"{name}\n")
            f.write(f"{value} {r} {g} {b} 255\n")


def import_dlabel(wb_command: str, ref_img: nib.Cifti2Image, values: np.ndarray, label_list: Path, out_path: Path) -> None:
    tmp = out_path.with_suffix("").with_suffix(".dscalar.nii")
    save_dscalar(ref_img, values, tmp)
    if not shutil.which(wb_command):
        print(f"[template_relabel] WARNING: {wb_command} unavailable; kept scalar {tmp}")
        return
    subprocess.run(
        [wb_command, "-cifti-label-import", str(tmp), str(label_list), str(out_path), "-discard-others"],
        check=True,
    )
    try:
        tmp.unlink()
    except OSError:
        pass


def dice_score(a: np.ndarray, b: np.ndarray) -> float:
    den = int(a.sum()) + int(b.sum())
    if den == 0:
        return 0.0
    return float(2 * int(np.logical_and(a, b).sum()) / den)


def relabel_density(
    communities: np.ndarray,
    cortex_idx: np.ndarray,
    cort_lookup: np.ndarray,
    spatial: np.ndarray,
    labels: Sequence[str],
    top_n: int,
    template_threshold: float,
    high_offset: int,
) -> Tuple[np.ndarray, List[dict]]:
    community_ids = np.unique(communities[communities > 0]).astype(int)
    n_net = len(labels)
    out = np.zeros(communities.size, dtype=np.float32)
    rows: List[dict] = []
    if community_ids.size == 0:
        return out, rows

    sizes = {int(cid): int(np.count_nonzero(communities == cid)) for cid in community_ids}
    ranked = sorted(community_ids.tolist(), key=lambda cid: (-sizes[int(cid)], int(cid)))
    eligible = set(int(cid) for cid in ranked[: max(0, int(top_n))])
    assigned_communities = set()
    assigned_templates = {}

    template_masks = spatial >= float(template_threshold)
    scores = {}
    for cid in community_ids:
        members = np.where(communities == cid)[0]
        cmask = cort_lookup[members] >= 0
        cort_members = members[cmask]
        cort_prior_rows = cort_lookup[cort_members]
        if cort_prior_rows.size:
            mean_prior = np.nanmean(spatial[cort_prior_rows, :], axis=0)
            comm_mask = np.zeros(spatial.shape[0], dtype=bool)
            comm_mask[cort_prior_rows] = True
            dice = np.array([dice_score(comm_mask, template_masks[:, k]) for k in range(n_net)], dtype=np.float64)
        else:
            mean_prior = np.zeros(n_net, dtype=np.float64)
            dice = np.zeros(n_net, dtype=np.float64)
        scores[int(cid)] = (mean_prior, dice)

    for template_idx in range(n_net):
        best = None
        for cid in ranked:
            cid = int(cid)
            if cid not in eligible or cid in assigned_communities:
                continue
            mean_prior, dice = scores[cid]
            key = (float(dice[template_idx]), float(mean_prior[template_idx]), sizes[cid] * -1, -cid)
            if best is None or key > best[0]:
                best = (key, cid)
        if best is not None:
            cid = best[1]
            assigned_communities.add(cid)
            assigned_templates[cid] = template_idx + 1

    for cid in ranked:
        cid = int(cid)
        new_id = int(assigned_templates.get(cid, high_offset + cid))
        out[communities == cid] = float(new_id)
        mean_prior, dice = scores[cid]
        template_id = int(assigned_templates.get(cid, 0))
        best_template = int(np.argmax(dice) + 1) if dice.size else 0
        rows.append(
            {
                "old_community_id": cid,
                "new_label_id": new_id,
                "community_size": sizes[cid],
                "rank_by_size": ranked.index(cid) + 1,
                "reserved_template_id": template_id,
                "reserved_template_name": labels[template_id - 1] if template_id else "",
                "best_template_id_by_dice": best_template,
                "best_template_name_by_dice": labels[best_template - 1] if best_template else "",
                "assigned_template_dice": float(dice[template_id - 1]) if template_id else 0.0,
                "assigned_template_mean_prior": float(mean_prior[template_id - 1]) if template_id else 0.0,
                "best_template_dice": float(np.max(dice)) if dice.size else 0.0,
                "best_template_mean_prior": float(mean_prior[best_template - 1]) if best_template else 0.0,
            }
        )
    return out, rows


def main() -> int:
    ap = argparse.ArgumentParser(description="Beta template-anchored relabeling for PFM Infomap communities")
    ap.add_argument("--communities-cifti", required=True)
    ap.add_argument("--priors-mat", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile-prefix", default="InfomapTemplateRelabelBeta")
    ap.add_argument("--top-n-communities", type=int, default=30)
    ap.add_argument("--template-threshold", type=float, default=0.05)
    ap.add_argument("--high-offset", type=int, default=1000)
    ap.add_argument("--wb-command", default="wb_command")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    img = nib.load(args.communities_cifti)
    data = np.asarray(img.get_fdata(), dtype=np.float64)
    if data.ndim != 2:
        raise ValueError("communities CIFTI must be 2D")

    brain_axis = img.header.get_axis(1)
    cortex_idx = cortical_grayordinates(brain_axis)
    cort_lookup = np.full(brain_axis.size, -1, dtype=np.int64)
    cort_lookup[cortex_idx] = np.arange(cortex_idx.size)

    spatial, labels, colors = load_priors(Path(args.priors_mat))
    if spatial.shape[0] < cortex_idx.size:
        raise ValueError(f"Priors.Spatial has {spatial.shape[0]} rows but CIFTI has {cortex_idx.size} cortical grayordinates")
    spatial = spatial[: cortex_idx.size, :]

    all_rows: List[dict] = []
    for density_idx in range(data.shape[0]):
        communities = np.rint(data[density_idx, :]).astype(int)
        relabeled, rows = relabel_density(
            communities,
            cortex_idx,
            cort_lookup,
            spatial,
            labels,
            int(args.top_n_communities),
            float(args.template_threshold),
            int(args.high_offset),
        )
        for row in rows:
            row["density_index"] = density_idx + 1
        dens_tag = f"Density{density_idx + 1:02d}"
        table = outdir / f"{args.outfile_prefix}_{dens_tag}_CommunityMap.csv"
        fields = [
            "density_index",
            "old_community_id",
            "new_label_id",
            "community_size",
            "rank_by_size",
            "reserved_template_id",
            "reserved_template_name",
            "best_template_id_by_dice",
            "best_template_name_by_dice",
            "assigned_template_dice",
            "assigned_template_mean_prior",
            "best_template_dice",
            "best_template_mean_prior",
        ]
        with table.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        label_list = outdir / f"{args.outfile_prefix}_{dens_tag}_LabelListFile.txt"
        write_label_list(label_list, rows, labels, colors)
        dlabel = outdir / f"{args.outfile_prefix}_{dens_tag}.dlabel.nii"
        import_dlabel(str(args.wb_command), img, relabeled, label_list, dlabel)
        all_rows.extend(rows)
        anchored = sum(1 for row in rows if int(row["reserved_template_id"]) > 0)
        print(f"[template_relabel] density={density_idx + 1} communities={len(rows)} anchored={anchored} wrote {dlabel}")

    summary = outdir / f"{args.outfile_prefix}_AllDensities_CommunityMap.csv"
    if all_rows:
        with summary.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)
    print(f"[template_relabel] wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Pure Python RidgePFM mapping for PFM."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.cifti2.cifti2_axes import SeriesAxis


def zscore_rows(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, keepdims=True)
    sd = np.where(sd > eps, sd, 1.0)
    return (x - mu) / sd


def bool_int(value: str) -> int:
    ivalue = int(value)
    if ivalue not in (0, 1):
        raise argparse.ArgumentTypeError("expected 0 or 1")
    return ivalue


def clamp_prob_smooth(p: np.ndarray, floor: float = 1e-6, alpha: float = 1e-5) -> np.ndarray:
    """Match the legacy MATLAB prior safety transform."""
    out = np.asarray(p, dtype=np.float32)
    if alpha > 0:
        out = (1.0 - float(alpha)) * out + float(alpha) / out.shape[1]
    return np.clip(out, float(floor), 1.0 - float(floor))


def spatial_prior_term(p: np.ndarray, transform: str) -> np.ndarray:
    if transform == "log":
        return np.log(np.maximum(p, 1e-6))
    if transform == "logit":
        q = clamp_prob_smooth(p)
        return np.log(q / (1.0 - q))
    raise ValueError(f"Unknown spatial prior transform: {transform}")


def save_scalar_like(ref_img: nib.Cifti2Image, values: np.ndarray, out_path: Path) -> None:
    out = values.reshape(-1, 1).astype(np.float32, copy=False)
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=1)
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(out.T, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def save_prob_dtseries(ref_img: nib.Cifti2Image, prob: np.ndarray, out_path: Path) -> None:
    # prob: n_gray x n_net -> dense timeseries laid out as maps x grayordinates
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=prob.shape[1])
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(prob.T.astype(np.float32), hdr, nifti_header=ref_img.nifti_header), str(out_path))


def write_label_list(out_txt: Path, labels: list[str], colors: np.ndarray) -> None:
    with out_txt.open("w") as f:
        for i, name in enumerate(labels, start=1):
            rgb = np.clip(np.rint(colors[i - 1] * 255.0), 0, 255).astype(int)
            f.write(f"{name}\n")
            f.write(f"{i} {rgb[0]} {rgb[1]} {rgb[2]} 255\n")


def parse_structures(axis, structures_csv: str) -> np.ndarray:
    requested = {s.strip().upper() for s in structures_csv.split(",") if s.strip()}
    keep = np.zeros((axis.size,), dtype=bool)
    for name, slc, _ in axis.iter_structures():
        stop = slc.stop if slc.stop is not None else axis.size
        tag = name.replace("CIFTI_STRUCTURE_", "").upper()
        if not requested or tag in requested:
            keep[slc.start:stop] = True
    return keep


def inject_subcortical_spatial_priors(
    axis,
    spatial_full: np.ndarray,
    subcort_priors_nii: str,
    n_net: int,
) -> None:
    if not subcort_priors_nii:
        return
    vol = nib.load(subcort_priors_nii)
    arr = vol.get_fdata(dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[..., np.newaxis]
    if arr.ndim != 4:
        raise ValueError(f"Expected 3D/4D NIfTI for subcortical priors: {subcort_priors_nii}, got {arr.shape}")
    if arr.shape[3] < n_net:
        raise ValueError(
            f"Subcortical priors maps ({arr.shape[3]}) must be >= network count ({n_net}): {subcort_priors_nii}"
        )

    used = 0
    for name, slc, bm in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            continue
        vox = np.asarray(bm.voxel, dtype=np.int64)
        if vox.size == 0:
            continue
        if (
            (vox[:, 0] < 0).any()
            or (vox[:, 1] < 0).any()
            or (vox[:, 2] < 0).any()
            or (vox[:, 0] >= arr.shape[0]).any()
            or (vox[:, 1] >= arr.shape[1]).any()
            or (vox[:, 2] >= arr.shape[2]).any()
        ):
            raise ValueError(f"Subcortical voxel indices are out of bounds for priors volume: {subcort_priors_nii}")
        spatial_full[slc, :] = arr[vox[:, 0], vox[:, 1], vox[:, 2], :n_net]
        used += int(vox.shape[0])

    print(f"[ridge] loaded ACPC subcortical priors: {subcort_priors_nii} (rows={used})")


def main() -> int:
    ap = argparse.ArgumentParser(description="Python RidgePFM")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--distance-npy", required=True, help="npy uint8 distance matrix")
    ap.add_argument("--priors-mat", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile", default="RidgePFM_VTX")
    ap.add_argument("--fc-weight", type=float, default=1.0)
    ap.add_argument("--fc-demean", type=bool_int, default=0, help="demean each target FC fingerprint after local edge exclusion")
    ap.add_argument("--positive-fc-only", type=bool_int, default=0, help="fit ridge model using only positive target FC edges")
    ap.add_argument("--spatial-weight", type=float, default=0.05)
    ap.add_argument(
        "--spatial-prior-transform",
        choices=("logit", "log"),
        default="logit",
        help="transform used for spatial prior fusion; logit matches legacy V4/V5 MATLAB RidgePFM",
    )
    ap.add_argument("--lambda", dest="lam", type=float, default=10.0)
    ap.add_argument("--local-exclusion-mm", type=float, default=10.0)
    ap.add_argument("--brain-structures-csv", default="")
    ap.add_argument("--subcort-priors-nii", default="", help="ACPC-space 4D subcortical priors NIfTI")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    img = nib.load(args.in_cifti)
    data_tg = np.asanyarray(img.dataobj).astype(np.float32, copy=False)  # time x gray
    data = data_tg.T  # gray x time
    axis = img.header.get_axis(1)

    # Cortex indices
    cortex_idx = []
    for name, slc, _ in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            stop = slc.stop if slc.stop is not None else axis.size
            cortex_idx.extend(range(slc.start, stop))
    cortex_idx = np.asarray(cortex_idx, dtype=np.int32)
    n_gray, n_time = data.shape
    n_cort = cortex_idx.size

    pri = sio.loadmat(args.priors_mat, squeeze_me=True, struct_as_record=False)["Priors"]
    pri_fc = np.asarray(pri.FC, dtype=np.float32)
    if pri_fc.shape[0] != n_cort:
        raise ValueError(f"Priors.FC rows ({pri_fc.shape[0]}) must match cortical count ({n_cort})")
    n_net = pri_fc.shape[1]
    pri_spatial = np.asarray(pri.Spatial, dtype=np.float32)
    if pri_spatial.shape[0] != n_cort:
        pri_spatial = pri_spatial[:n_cort, :]
    labels_raw = np.asarray(pri.NetworkLabels).ravel()
    labels = [str(x).strip() for x in labels_raw[:n_net]]
    colors = np.asarray(pri.NetworkColors, dtype=np.float32)[:n_net, :]

    # Selected grayordinates by structure
    selected = parse_structures(axis, args.brain_structures_csv)
    if not selected.any():
        selected[:] = True
    good_idx = np.where(selected)[0]

    # Prepare design
    x_cort = zscore_rows(data[cortex_idx, :].astype(np.float64))
    a = pri_fc.astype(np.float64)
    a = a - a.mean(axis=0, keepdims=True)
    a /= np.maximum(np.sqrt((a * a).sum(axis=0, keepdims=True)), 1e-8)
    w = np.linalg.solve(a.T @ a + args.lam * np.eye(n_net), a.T)  # n_net x n_cort

    # Distance matrix (for local exclusion on cortical targets only)
    dist = np.load(args.distance_npy, mmap_mode="r")
    d_cort = dist[np.ix_(cortex_idx, cortex_idx)].astype(np.float32)

    # Spatial priors full
    spatial_full = np.full((n_gray, n_net), 0.5, dtype=np.float32)
    spatial_full[cortex_idx, :] = pri_spatial
    inject_subcortical_spatial_priors(axis, spatial_full, args.subcort_priors_nii, n_net)
    spatial_full = clamp_prob_smooth(spatial_full)

    label_idx = np.zeros((n_gray,), dtype=np.int32)
    prob_best = np.zeros((n_gray,), dtype=np.float32)
    r2 = np.zeros((n_gray,), dtype=np.float32)
    prob_all = np.zeros((n_gray, n_net), dtype=np.float32)

    for i in range(0, good_idx.size, args.chunk):
        j = min(i + args.chunk, good_idx.size)
        g = good_idx[i:j]
        y = zscore_rows(data[g, :].astype(np.float64))
        m = (x_cort @ y.T) / max(n_time - 1, 1)  # n_cort x chunk

        # Demean FC fingerprints after local exclusion; useful without GSR/MGTR.
        for c, gi in enumerate(g):
            keep = np.ones((n_cort,), dtype=bool)
            pos = np.where(cortex_idx == gi)[0]
            if pos.size == 1:
                keep &= d_cort[:, pos[0]] > args.local_exclusion_mm
            if int(args.fc_demean) == 1:
                if keep.any():
                    m[:, c] -= float(m[keep, c].mean())
                else:
                    m[:, c] = 0.0
            m[~keep, c] = 0.0

        if int(args.positive_fc_only) == 1:
            beta = np.zeros((g.size, n_net), dtype=np.float64)
            for c in range(g.size):
                keep_pos = m[:, c] > 0
                if not keep_pos.any():
                    continue
                ak = a[keep_pos, :]
                beta[c, :] = np.linalg.solve(ak.T @ ak + args.lam * np.eye(n_net), ak.T @ m[keep_pos, c])
        else:
            beta = (w @ m).T  # chunk x n_net
        beta_mu = beta.mean(axis=1, keepdims=True)
        beta_sd = np.maximum(beta.std(axis=1, keepdims=True), 1e-8)
        beta_z = (beta - beta_mu) / beta_sd

        score = args.fc_weight * beta_z + args.spatial_weight * spatial_prior_term(
            spatial_full[g, :], args.spatial_prior_transform
        )
        score -= score.max(axis=1, keepdims=True)
        p = np.exp(score)
        p /= np.maximum(p.sum(axis=1, keepdims=True), 1e-8)

        best = np.argmax(p, axis=1)
        label_idx[g] = best + 1
        prob_best[g] = p[np.arange(p.shape[0]), best]
        prob_all[g, :] = p

        # crude R2 on FC map fit
        m_hat = (a @ beta.T)  # n_cort x chunk
        ss_res = ((m - m_hat) ** 2).sum(axis=0)
        m_mu = m.mean(axis=0, keepdims=True)
        ss_tot = ((m - m_mu) ** 2).sum(axis=0) + 1e-8
        r2[g] = (1.0 - (ss_res / ss_tot)).astype(np.float32)

        if (i // args.chunk + 1) % 10 == 0:
            print(f"[ridge] processed {j}/{good_idx.size}")

    # Write outputs
    label_tmp = outdir / "Tmp_labels.dtseries.nii"
    save_scalar_like(img, label_idx.astype(np.float32), label_tmp)
    labfile = outdir / "LabelListFile.txt"
    write_label_list(labfile, labels, colors)
    dlabel = outdir / f"{args.outfile}.dlabel.nii"
    subprocess.run(
        ["wb_command", "-cifti-label-import", str(label_tmp), str(labfile), str(dlabel), "-discard-others"],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(dlabel), "-border", args.left_surf, str(outdir / f"{args.outfile}.L.border")],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(dlabel), "-border", args.right_surf, str(outdir / f"{args.outfile}.R.border")],
        check=True,
    )

    save_scalar_like(img, r2, outdir / f"{args.outfile}_R2.dtseries.nii")
    save_prob_dtseries(img, prob_all, outdir / f"{args.outfile}_ProbMaps.dtseries.nii")
    try:
        if label_tmp.exists():
            label_tmp.unlink()
        if labfile.exists():
            labfile.unlink()
    except Exception:
        pass
    print(f"[ridge] wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

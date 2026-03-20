#!/usr/bin/env python3
"""Simple Python areal parcellation pass from ridge WTA labels.

This implementation creates per-network connected components over full
grayordinate neighbors and writes a parcel dlabel output.
"""

import argparse
import subprocess
from collections import deque
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.affines import apply_affine
from nibabel.cifti2.cifti2_axes import SeriesAxis


def save_scalar_like(ref_img: nib.Cifti2Image, values: np.ndarray, out_path: Path) -> None:
    out = values.reshape(-1, 1).astype(np.float32, copy=False)
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=1)
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(out.T, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def build_full_neighbors(img: nib.Cifti2Image, neighbors_cortex: np.ndarray) -> np.ndarray:
    bm = img.header.get_axis(1)
    n = bm.size
    coords = np.full((n, 3), np.nan, dtype=np.float64)
    cortex_mask = np.zeros((n,), dtype=bool)
    left_right_idx = []

    for name, slc, sub in bm.iter_structures():
        stop = slc.stop if slc.stop is not None else n
        gidx = np.arange(slc.start, stop, dtype=np.int32)
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            cortex_mask[gidx] = True
            left_right_idx.extend(gidx.tolist())
        vox = np.asarray(sub.voxel, dtype=np.float64)
        if vox.shape[0] == gidx.size:
            coords[gidx, :] = apply_affine(bm.affine, vox)

    idx_c = np.asarray(left_right_idx, dtype=np.int32)
    idx_s = np.setdiff1d(np.arange(n, dtype=np.int32), idx_c, assume_unique=True)
    k = max(neighbors_cortex.shape[1] + 1, 7)
    out = np.zeros((n, k), dtype=np.int32)
    out[idx_c, 0] = idx_c

    surf = neighbors_cortex
    if surf.shape[0] == idx_c.size and np.all(surf[:, 0] == np.arange(1, surf.shape[0] + 1)):
        surf = surf[:, 1:]
    mask = surf > 0
    surf_m = surf.copy()
    surf_m[mask] = idx_c[surf[mask] - 1]
    cpy = min(surf_m.shape[1], k - 1)
    out[idx_c, 1 : 1 + cpy] = surf_m[:, :cpy]

    # Subcortical 6-neighborhood from coords (axis-aligned unit steps inferred)
    spos = coords[idx_s, :]
    valid = np.all(np.isfinite(spos), axis=1)
    steps = []
    for d in range(3):
        vals = np.unique(spos[valid, d])
        steps.append(np.min(np.diff(vals)) if vals.size >= 2 else 1.0)
    dirs = np.array(
        [
            [steps[0], 0, 0],
            [-steps[0], 0, 0],
            [0, steps[1], 0],
            [0, -steps[1], 0],
            [0, 0, steps[2]],
            [0, 0, -steps[2]],
        ],
        dtype=np.float64,
    )
    lookup = {}
    for i, gi in enumerate(idx_s):
        if valid[i]:
            lookup[tuple(np.round(spos[i, :], 6))] = gi
    out[idx_s, 0] = idx_s
    for i, gi in enumerate(idx_s):
        if not valid[i]:
            continue
        base = spos[i, :]
        nbrs = []
        for d in dirs:
            key = tuple(np.round(base + d, 6))
            if key in lookup:
                nbrs.append(lookup[key])
        if nbrs:
            take = min(len(nbrs), k - 1)
            out[gi, 1 : 1 + take] = np.asarray(nbrs[:take], dtype=np.int32)
    return out


def connected_components(labels: np.ndarray, neighbors: np.ndarray, min_size: int) -> np.ndarray:
    n = labels.size
    visited = np.zeros((n,), dtype=bool)
    parcels = np.zeros((n,), dtype=np.int32)
    next_id = 1
    for net in np.unique(labels):
        if net <= 0:
            continue
        idx = np.where(labels == net)[0]
        idx_set = set(idx.tolist())
        for s in idx:
            if visited[s]:
                continue
            q = deque([int(s)])
            comp = []
            visited[s] = True
            while q:
                u = q.popleft()
                comp.append(u)
                for v in neighbors[u, 1:]:
                    v = int(v)
                    if v <= 0 or v >= n:
                        continue
                    if v in idx_set and not visited[v]:
                        visited[v] = True
                        q.append(v)
            if len(comp) >= min_size:
                parcels[np.asarray(comp, dtype=np.int32)] = next_id
                next_id += 1
    return parcels


def main() -> int:
    ap = argparse.ArgumentParser(description="Simple areal parcellation from ridge WTA labels")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--wta-dlabel", required=True)
    ap.add_argument("--neighbors-mat", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile", default="RidgeFusion_VTX+ArealParcellation")
    ap.add_argument("--min-size", type=int, default=10)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cimg = nib.load(args.in_cifti)
    limg = nib.load(args.wta_dlabel)
    labels = np.asanyarray(limg.dataobj).squeeze().astype(np.int32, copy=False)

    nmat = sio.loadmat(args.neighbors_mat, squeeze_me=True, struct_as_record=False)
    nkey = [k for k in nmat.keys() if not k.startswith("__")][0]
    neighbors_cortex = np.asarray(nmat[nkey], dtype=np.int32)
    full_neighbors = build_full_neighbors(cimg, neighbors_cortex)
    parcels = connected_components(labels, full_neighbors, min_size=args.min_size)

    tmp = outdir / "Tmp_Areal.dtseries.nii"
    save_scalar_like(cimg, parcels.astype(np.float32), tmp)
    out_dlabel = outdir / f"{args.outfile}.dlabel.nii"
    # Reuse simple random-ish palette
    uniq = np.unique(parcels)
    labfile = outdir / "LabelListFile_Areal.txt"
    with labfile.open("w") as f:
        for u in uniq:
            if u <= 0:
                continue
            f.write(f"Parcel_{int(u)}\n")
            r = (37 * int(u)) % 255
            g = (73 * int(u)) % 255
            b = (109 * int(u)) % 255
            f.write(f"{int(u)} {r} {g} {b} 255\n")
    subprocess.run(
        ["wb_command", "-cifti-label-import", str(tmp), str(labfile), str(out_dlabel), "-discard-others"],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(out_dlabel), "-border", args.left_surf, str(outdir / f"{args.outfile}.L.border")],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(out_dlabel), "-border", args.right_surf, str(outdir / f"{args.outfile}.R.border")],
        check=True,
    )
    if tmp.exists():
        tmp.unlink()
    if labfile.exists():
        labfile.unlink()
    print(f"[areal] wrote {out_dlabel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

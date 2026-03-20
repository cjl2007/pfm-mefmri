#!/usr/bin/env python3
"""Build PFM distance matrix (pure Python + wb_command).

Distance model (single default):
- Cortex within hemisphere: geodesic distance on midthickness surface
- Cortex across hemispheres: Euclidean distance in 3D
- Subcortex to all grayordinates: Euclidean distance in 3D

Output is uint8 `.npy` matrix (n_gray x n_gray), rounded and clipped to [0,255].
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.affines import apply_affine
from numpy.lib.format import open_memmap


def to_uint8_dist(x: np.ndarray) -> np.ndarray:
    return np.clip(np.rint(x), 0, 255).astype(np.uint8, copy=False)


def load_gray_coords(ref_cifti: Path, left_surf: Path, right_surf: Path):
    img = nib.load(str(ref_cifti))
    bm = img.header.get_axis(1)
    n_gray = bm.size

    l_coords = np.asarray(nib.load(str(left_surf)).darrays[0].data, dtype=np.float64)
    r_coords = np.asarray(nib.load(str(right_surf)).darrays[0].data, dtype=np.float64)

    coords = np.full((n_gray, 3), np.nan, dtype=np.float64)
    left_idx, right_idx = [], []
    left_vert, right_vert = [], []

    for name, slc, sub in bm.iter_structures():
        stop = slc.stop if slc.stop is not None else n_gray
        gidx = np.arange(slc.start, stop, dtype=np.int32)
        if name == "CIFTI_STRUCTURE_CORTEX_LEFT":
            v = np.asarray(sub.vertex, dtype=np.int32)
            coords[gidx, :] = l_coords[v, :]
            left_idx = gidx.tolist()
            left_vert = v.tolist()
        elif name == "CIFTI_STRUCTURE_CORTEX_RIGHT":
            v = np.asarray(sub.vertex, dtype=np.int32)
            coords[gidx, :] = r_coords[v, :]
            right_idx = gidx.tolist()
            right_vert = v.tolist()
        else:
            vox = np.asarray(sub.voxel, dtype=np.float64)
            coords[gidx, :] = apply_affine(bm.affine, vox)

    if np.isnan(coords).any():
        raise RuntimeError("Failed to derive full grayordinate coordinate set")

    left_idx = np.asarray(left_idx, dtype=np.int32)
    right_idx = np.asarray(right_idx, dtype=np.int32)
    left_vert = np.asarray(left_vert, dtype=np.int32)
    right_vert = np.asarray(right_vert, dtype=np.int32)
    sub_idx = np.setdiff1d(np.arange(n_gray, dtype=np.int32), np.r_[left_idx, right_idx], assume_unique=True)
    return n_gray, coords, left_idx, right_idx, left_vert, right_vert, sub_idx


def geodesic_all_to_all(surf: Path, out_dconn: Path) -> None:
    subprocess.run(
        ["wb_command", "-surface-geodesic-distance-all-to-all", str(surf), str(out_dconn), "-limit", "255"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def dconn_to_uint8_submatrix(dconn_path: Path, row_verts: np.ndarray, col_verts: np.ndarray) -> np.ndarray:
    img = nib.load(str(dconn_path))
    data = np.asanyarray(img.dataobj)  # may be memmapped
    sub = data[np.ix_(row_verts, col_verts)].astype(np.float64, copy=False)
    sub[sub < 0] = 255.0  # wb -limit emits -1 for beyond limit
    return to_uint8_dist(sub)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build PFM distance matrix with geodesic+euclidean defaults")
    ap.add_argument("--ref-cifti", required=True)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    ap.add_argument("--out-npy", required=True)
    ap.add_argument("--chunk-rows", type=int, default=128)
    args = ap.parse_args()

    ref_cifti = Path(args.ref_cifti)
    left_surf = Path(args.left_surf)
    right_surf = Path(args.right_surf)
    out_npy = Path(args.out_npy)
    out_npy.parent.mkdir(parents=True, exist_ok=True)

    n_gray, coords, left_idx, right_idx, left_vert, right_vert, sub_idx = load_gray_coords(
        ref_cifti, left_surf, right_surf
    )
    print(f"[dist] n_gray={n_gray}, n_left={left_idx.size}, n_right={right_idx.size}, n_sub={sub_idx.size}")

    D = open_memmap(str(out_npy), mode="w+", dtype=np.uint8, shape=(n_gray, n_gray))
    D[:] = 255

    with tempfile.TemporaryDirectory(prefix="pfm_geo_all2all_") as tdir:
        tdir = Path(tdir)
        l_dconn = tdir / "left_geodesic.dconn.nii"
        r_dconn = tdir / "right_geodesic.dconn.nii"
        print("[dist] computing left hemisphere all-to-all geodesic")
        geodesic_all_to_all(left_surf, l_dconn)
        print("[dist] computing right hemisphere all-to-all geodesic")
        geodesic_all_to_all(right_surf, r_dconn)

        print("[dist] extracting cortex within-hemisphere blocks")
        D[np.ix_(left_idx, left_idx)] = dconn_to_uint8_submatrix(l_dconn, left_vert, left_vert)
        D[np.ix_(right_idx, right_idx)] = dconn_to_uint8_submatrix(r_dconn, right_vert, right_vert)

    print("[dist] filling cortex cross-hemisphere Euclidean block")
    chunk = int(args.chunk_rows)
    lxyz = coords[left_idx]
    rxyz = coords[right_idx]
    for i in range(0, left_idx.size, chunk):
        j = min(i + chunk, left_idx.size)
        d = np.linalg.norm(lxyz[i:j, None, :] - rxyz[None, :, :], axis=2)
        du = to_uint8_dist(d)
        D[np.ix_(left_idx[i:j], right_idx)] = du
        D[np.ix_(right_idx, left_idx[i:j])] = du.T

    print("[dist] filling subcortical Euclidean rows/cols")
    allxyz = coords
    for i in range(0, sub_idx.size, chunk):
        j = min(i + chunk, sub_idx.size)
        rows = sub_idx[i:j]
        d = np.linalg.norm(allxyz[rows, None, :] - allxyz[None, :, :], axis=2)
        du = to_uint8_dist(d)
        D[rows, :] = du
        D[:, rows] = du.T
        if (i // chunk + 1) % 20 == 0:
            print(f"[dist] sub rows done: {j}/{sub_idx.size}")

    # enforce exact symmetry + zero diagonal
    print("[dist] finalizing symmetry")
    for i in range(0, n_gray, chunk):
        j = min(i + chunk, n_gray)
        blk = D[i:j, :]
        blk_t = D[:, i:j].T
        D[i:j, :] = np.minimum(blk, blk_t)
    np.fill_diagonal(D, 0)
    D.flush()
    print(f"[dist] wrote {out_npy}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


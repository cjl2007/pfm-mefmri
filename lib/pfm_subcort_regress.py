#!/usr/bin/env python3
"""Regress nearby cortical signal from subcortical voxel timeseries.

Python port of the legacy regress_cortical_signals.m behavior:
  - Identify subcortical voxels adjacent to cortex within a distance threshold
  - For each adjacent subcortical voxel: regress y ~ mean(nearby cortex) + intercept
  - Replace only regressed subcortical rows with residuals
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import nibabel as nib
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from nibabel.affines import apply_affine
from nibabel.cifti2.cifti2_axes import BrainModelAxis, SeriesAxis


def load_coords_and_masks(
    img: nib.Cifti2Image, left_surf: Path, right_surf: Path
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    bm = img.header.get_axis(1)
    if not isinstance(bm, BrainModelAxis):
        raise ValueError("Expected BrainModelAxis as axis 1")

    l_surf = nib.load(str(left_surf))
    r_surf = nib.load(str(right_surf))
    l_coords = np.asarray(l_surf.darrays[0].data, dtype=np.float64)
    r_coords = np.asarray(r_surf.darrays[0].data, dtype=np.float64)

    n_gray = bm.size
    coords = np.full((n_gray, 3), np.nan, dtype=np.float64)
    cortex = np.zeros((n_gray,), dtype=bool)

    for name, slc, sub in bm.iter_structures():
        stop = slc.stop if slc.stop is not None else n_gray
        gidx = np.arange(slc.start, stop, dtype=np.int32)
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            cortex[gidx] = True
            vert = np.asarray(sub.vertex, dtype=np.int32)
            if name == "CIFTI_STRUCTURE_CORTEX_LEFT":
                coords[gidx, :] = l_coords[vert, :]
            else:
                coords[gidx, :] = r_coords[vert, :]
            continue
        vox = np.asarray(sub.voxel, dtype=np.float64)
        coords[gidx, :] = apply_affine(bm.affine, vox)

    if np.isnan(coords).any():
        raise RuntimeError("Could not derive coordinates for all grayordinates")

    subcort = ~cortex
    return coords, cortex, subcort


def save_dtseries(data_time_by_gray: np.ndarray, ref_img: nib.Cifti2Image, out_path: Path) -> None:
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    if len(axes) != 2 or not isinstance(axes[0], SeriesAxis):
        raise ValueError("Expected dtseries with SeriesAxis first")
    new_series = SeriesAxis(axes[0].start, axes[0].step, data_time_by_gray.shape[0], unit=axes[0].unit)
    new_header = nib.Cifti2Header.from_axes((new_series, axes[1]))
    out_img = nib.Cifti2Image(data_time_by_gray.astype(np.float32), new_header, nifti_header=ref_img.nifti_header)
    nib.save(out_img, str(out_path))


def compute_regression_masks_from_distance(
    row: int,
    distance_row: np.ndarray,
    cortex: np.ndarray,
    distance_mm: float,
) -> tuple[bool, np.ndarray]:
    cortex_idx = np.where(cortex)[0]
    near_cortex = np.asarray(distance_row[cortex_idx] <= distance_mm)
    should_regress = bool(near_cortex.any())
    nb_mask = np.zeros_like(cortex, dtype=bool)
    nb_mask[cortex_idx] = near_cortex
    return should_regress, nb_mask


def regress_runs(y, neighbor_mean, runs):
    """Separate run slopes/intercepts; inputs are the retained fitting frames."""
    design = np.column_stack((neighbor_mean, np.ones(len(y))))
    residual = np.asarray(y, dtype=np.float64).copy()
    for rows in runs:
        if len(rows) < 3:
            raise ValueError('Each run needs at least three retained frames for regression')
        beta = np.linalg.lstsq(design[rows], y[rows], rcond=None)[0]
        residual[rows] = y[rows] - design[rows] @ beta
    return residual


def main() -> int:
    ap = argparse.ArgumentParser(description="Subcortical neighborhood regression for CIFTI dtseries")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--out-cifti", required=True)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    ap.add_argument("--distance-npy", default="", help="Optional checked full grayordinate distance matrix; default derives physical neighborhoods from surfaces and CIFTI voxels")
    ap.add_argument("--distance-mm", type=float, default=20.0)
    ap.add_argument("--scan-idx", default="", help="Frame-aligned MAT ScanIdx [session, run]; fit separate run intercepts/slopes")
    ap.add_argument("--max-subvox", type=int, default=0, help="Optional debug cap on #subcortical voxels (0=all)")
    args = ap.parse_args()

    in_path = Path(args.in_cifti)
    out_path = Path(args.out_cifti)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(in_path))
    data = np.asanyarray(img.dataobj).astype(np.float32, copy=False)  # time x grayordinates
    if data.ndim != 2:
        raise ValueError(f"Expected 2D dtseries (time x grayordinates), got {data.shape}")

    coords, cortex, subcort = load_coords_and_masks(
        img, left_surf=Path(args.left_surf), right_surf=Path(args.right_surf)
    )
    sub_idx = np.where(subcort)[0]
    if args.max_subvox > 0:
        sub_idx = sub_idx[: args.max_subvox]

    out = data.copy()
    if args.scan_idx:
        groups = np.asarray(loadmat(args.scan_idx)['ScanIdx'])
        if groups.shape != (len(data), 2) or not np.isfinite(groups).all():
            raise ValueError('ScanIdx must align exactly with input frames')
        runs = [np.flatnonzero((groups == pair).all(1)) for pair in np.unique(groups, axis=0)]
    else:
        runs = [np.arange(len(data))]
    dist_mm = float(args.distance_mm)
    dist = np.load(args.distance_npy, mmap_mode="r") if args.distance_npy else None
    if dist is not None and dist.shape != (data.shape[1], data.shape[1]):
        raise ValueError(f"Distance matrix shape {dist.shape} does not match CIFTI grayordinates {data.shape[1]}")

    n_regressed = 0
    n_skipped = 0
    cortex_idx = np.flatnonzero(cortex)
    tree = cKDTree(coords[cortex_idx])
    for i, gi in enumerate(sub_idx, start=1):
        if dist is not None:
            should_regress, nb_mask = compute_regression_masks_from_distance(
                int(gi),
                np.asarray(dist[gi, :]),
                cortex,
                dist_mm,
            )
        else:
            # Preserve rounded-distance <= threshold convention without the
            # legacy one-voxel shift or dependence on an old distance matrix.
            near = np.asarray(tree.query_ball_point(coords[gi], dist_mm + 0.5), dtype=int)
            near = near[np.linalg.norm(coords[cortex_idx[near]] - coords[gi], axis=1) < dist_mm + 0.5]
            nb_mask = np.zeros_like(cortex)
            nb_mask[cortex_idx[near]] = True
            should_regress = bool(len(near))

        if not should_regress or not nb_mask.any():
            n_skipped += 1
            continue

        nb_ts = data[:, nb_mask]
        if nb_ts.shape[1] > 1:
            nb_mean = nb_ts.mean(axis=1, dtype=np.float64)
        else:
            nb_mean = nb_ts[:, 0].astype(np.float64, copy=False)

        y = data[:, gi].astype(np.float64, copy=False)
        resid = regress_runs(y, nb_mean, runs)
        out[:, gi] = resid.astype(np.float32, copy=False)
        n_regressed += 1

        if i % 1000 == 0:
            print(f"[subcort_regress] processed {i}/{len(sub_idx)} subcortical voxels")

    save_dtseries(out, img, out_path)
    print(f"[subcort_regress] regressed={n_regressed} skipped={n_skipped} neighbor_source=cortex")
    print(f"[subcort_regress] wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

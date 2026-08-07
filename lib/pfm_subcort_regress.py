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


def main() -> int:
    ap = argparse.ArgumentParser(description="Subcortical neighborhood regression for CIFTI dtseries")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--out-cifti", required=True)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    ap.add_argument("--distance-npy", default="", help="Full grayordinate distance matrix; preferred for legacy-compatible neighborhoods")
    ap.add_argument("--distance-mm", type=float, default=20.0)
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
    ones = np.ones((data.shape[0], 1), dtype=np.float64)
    dist_mm = float(args.distance_mm)
    dist = np.load(args.distance_npy, mmap_mode="r") if args.distance_npy else None
    if dist is not None and dist.shape != (data.shape[1], data.shape[1]):
        raise ValueError(f"Distance matrix shape {dist.shape} does not match CIFTI grayordinates {data.shape[1]}")

    n_regressed = 0
    n_skipped = 0
    for i, gi in enumerate(sub_idx, start=1):
        if dist is not None:
            should_regress, nb_mask = compute_regression_masks_from_distance(
                int(gi),
                np.asarray(dist[gi, :]),
                cortex,
                dist_mm,
            )
        else:
            d = np.linalg.norm(coords - coords[gi][None, :], axis=1)
            should_regress = bool(np.any(d[cortex] <= dist_mm))
            nb_mask = cortex & (d <= dist_mm)

        if not should_regress or not nb_mask.any():
            n_skipped += 1
            continue

        nb_ts = data[:, nb_mask]
        if nb_ts.shape[1] > 1:
            nb_mean = nb_ts.mean(axis=1, dtype=np.float64)
        else:
            nb_mean = nb_ts[:, 0].astype(np.float64, copy=False)

        y = data[:, gi].astype(np.float64, copy=False)
        X = np.column_stack((nb_mean, ones[:, 0]))
        beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
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

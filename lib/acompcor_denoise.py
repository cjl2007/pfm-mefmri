#!/usr/bin/env python3
"""Standalone beta aCompCor denoising module.

Implements a Python-only nuisance-regression workflow inspired by the
description provided by the user:
  - WM and ventricular nuisance masks from FreeSurfer aseg, resampled to fMRI.
  - Extra-axial nuisance mask from temporal SD% outside a dilated brain mask.
  - Adaptive per-compartment dimensionality reduction using an eigenvalue
    condition-number threshold.
  - Second-stage SVD conditioning on the combined nuisance design matrix.
  - Regression fit using only uncensored frames; censored frames are excised
    from the main output after denoising.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy.ndimage import binary_dilation, binary_erosion


WM_LABELS = np.array([2, 7, 41, 46, 77, 78, 79, 251, 252, 253, 254, 255], dtype=np.int32)
VENT_LABELS = np.array([4, 5, 14, 15, 24, 31, 43, 44, 63], dtype=np.int32)


@dataclass
class CompartmentResult:
    name: str
    mask: np.ndarray
    n_voxels: int
    n_components: int
    eigenvalues: List[float]
    regressors: np.ndarray


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Run beta aCompCor denoising on a 4D fMRI series.")
    ap.add_argument("--subdir", required=True, help="Subject directory")
    ap.add_argument("--input", required=True, help="Input 4D fMRI timeseries")
    ap.add_argument("--motion", required=True, help="6-column motion parameter file")
    ap.add_argument("--output-base", required=True, help="Output base path for censored-cleaned NIfTI")
    ap.add_argument("--output-full-base", default="", help="Optional output base path for full-length cleaned NIfTI")
    ap.add_argument("--censor", default="", help="Optional text file containing frame keep/censor flags")
    ap.add_argument("--censor-keep-value", type=int, default=1, help="Value in censor file that means keep frame (default: 1)")
    ap.add_argument("--brain-mask", default="", help="Optional explicit functional-space brain mask")
    ap.add_argument("--ribbon-mask", default="", help="Optional explicit functional-space cortical ribbon mask")
    ap.add_argument("--aseg", default="", help="Optional explicit FreeSurfer aseg/aparc+aseg volume")
    ap.add_argument("--extraaxial-std-pct", type=float, default=2.5, help="Temporal SD percentage threshold for extra-axial mask")
    ap.add_argument("--extraaxial-dilate-iters", type=int, default=1, help="Brain-mask dilation iterations before extra-axial exclusion")
    ap.add_argument("--wm-erode-iters", type=int, default=2, help="WM mask erosion iterations after resampling")
    ap.add_argument("--compartment-cond-max", type=float, default=30.0, help="Max condition number per nuisance compartment")
    ap.add_argument("--design-cond-max", type=float, default=250.0, help="Max condition number for combined nuisance design")
    return ap.parse_args()


def load_nifti(path: Path) -> nib.Nifti1Image:
    if not path.is_file():
        raise FileNotFoundError(f"Missing NIfTI: {path}")
    return nib.load(str(path))


def choose_default_paths(subdir: Path, args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    xfms_dir = os.environ.get("FUNC_XFMS_DIRNAME", "rest")
    atlas_space = os.environ.get("AtlasSpace", "T1w")

    if args.brain_mask:
        brain_mask = Path(args.brain_mask)
    else:
        brain_mask = subdir / "func" / "xfms" / xfms_dir / "T1w_acpc_brain_func_mask.nii.gz"
        if not brain_mask.is_file():
            brain_mask = subdir / "func" / "xfms" / xfms_dir / "T1w_func_brain_mask.nii.gz"

    if args.ribbon_mask:
        ribbon_mask = Path(args.ribbon_mask)
    else:
        ribbon_name = "CorticalRibbon_nonlin_func_mask.nii.gz" if atlas_space == "MNINonlinear" else "CorticalRibbon_acpc_func_mask.nii.gz"
        ribbon_mask = subdir / "func" / "xfms" / xfms_dir / ribbon_name

    if args.aseg:
        aseg = Path(args.aseg)
    else:
        aseg = subdir / "anat" / "T1w" / "aparc+aseg.nii.gz"

    return brain_mask, ribbon_mask, aseg


def load_data_matrix(img: nib.Nifti1Image) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    data = np.asanyarray(img.dataobj, dtype=np.float64)
    if data.ndim != 4:
        raise RuntimeError(f"Expected 4D input, got shape {data.shape}")
    shape = tuple(int(v) for v in data.shape)  # type: ignore[assignment]
    return data.reshape((-1, shape[3])), shape


def load_motion(path: Path, n_t: int) -> np.ndarray:
    motion = np.loadtxt(path, dtype=np.float64)
    if motion.ndim == 1:
        motion = motion[:, None]
    if motion.shape[0] != n_t and motion.shape[1] == n_t:
        motion = motion.T
    if motion.shape[0] != n_t:
        raise RuntimeError(f"Motion timepoints mismatch: got {motion.shape[0]}, expected {n_t}")
    if motion.shape[1] != 6:
        raise RuntimeError(f"Expected 6-column motion file, got shape {motion.shape}: {path}")
    return motion


def load_keep_mask(path: Optional[Path], n_t: int, keep_value: int) -> np.ndarray:
    if path is None:
        return np.ones(n_t, dtype=bool)
    vals = np.loadtxt(path, dtype=np.float64)
    vals = np.ravel(vals)
    if vals.shape[0] != n_t:
        raise RuntimeError(f"Censor vector length mismatch: got {vals.shape[0]}, expected {n_t}")
    return vals == keep_value


def resample_label_mask(aseg_img: nib.Nifti1Image, ref_img: nib.Nifti1Image, labels: np.ndarray) -> np.ndarray:
    # Resample against spatial axes only. Passing a 4D target image to
    # nibabel can trigger a 5x5 affine path in older nibabel versions.
    to_vox_map = (tuple(int(v) for v in ref_img.shape[:3]), ref_img.affine)
    aseg_res = resample_from_to(aseg_img, to_vox_map, order=0)
    aseg = np.asanyarray(aseg_res.dataobj)
    return np.isin(np.rint(aseg).astype(np.int32), labels)


def clean_binary_mask(mask: np.ndarray, erode_iters: int = 0, dilate_iters: int = 0) -> np.ndarray:
    out = mask.astype(bool)
    if erode_iters > 0:
        out = binary_erosion(out, iterations=erode_iters)
    if dilate_iters > 0:
        out = binary_dilation(out, iterations=dilate_iters)
    return out


def safe_standardize_columns(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0:
        return x
    mu = np.nanmean(x, axis=0, keepdims=True)
    sd = np.nanstd(x, axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    out = (x - mu) / sd
    out[~np.isfinite(out)] = 0.0
    return out


def compute_global_signal(data_2d: np.ndarray, ribbon_mask_3d: np.ndarray) -> np.ndarray:
    idx = np.flatnonzero(ribbon_mask_3d.reshape(-1) > 0)
    if idx.size == 0:
        raise RuntimeError("Cortical ribbon mask is empty.")
    return np.nanmean(data_2d[idx, :], axis=0)


def temporal_sd_percent(data_2d: np.ndarray) -> np.ndarray:
    mean = np.nanmean(data_2d, axis=1)
    sd = np.nanstd(data_2d, axis=1, ddof=1)
    denom = np.where(np.abs(mean) > 1e-8, np.abs(mean), np.nan)
    sd_pct = 100.0 * sd / denom
    sd_pct[~np.isfinite(sd_pct)] = 0.0
    return sd_pct


def compartment_regressors(
    data_2d: np.ndarray,
    mask_3d: np.ndarray,
    cond_max: float,
    name: str,
) -> CompartmentResult:
    idx = np.flatnonzero(mask_3d.reshape(-1) > 0)
    if idx.size == 0:
        return CompartmentResult(name=name, mask=mask_3d, n_voxels=0, n_components=0, eigenvalues=[], regressors=np.zeros((data_2d.shape[1], 0)))

    y = data_2d[idx, :].T
    y = y - np.nanmean(y, axis=0, keepdims=True)
    y[:, np.nanstd(y, axis=0) == 0] = 0.0

    u, s, _ = np.linalg.svd(y, full_matrices=False)
    eigvals = (s ** 2).tolist()
    positive = s > 0
    if not np.any(positive):
        return CompartmentResult(name=name, mask=mask_3d, n_voxels=idx.size, n_components=0, eigenvalues=eigvals, regressors=np.zeros((data_2d.shape[1], 0)))

    keep = positive & ((s[0] / np.maximum(s, 1e-12)) <= np.sqrt(cond_max))
    if not np.any(keep):
        keep[0] = True
    regressors = u[:, keep]
    return CompartmentResult(name=name, mask=mask_3d, n_voxels=idx.size, n_components=int(regressors.shape[1]), eigenvalues=eigvals, regressors=regressors)


def condition_design(x: np.ndarray, cond_max: float) -> Tuple[np.ndarray, List[float], int]:
    if x.size == 0:
        return x, [], 0
    xz = safe_standardize_columns(x)
    u, s, _ = np.linalg.svd(xz, full_matrices=False)
    eigvals = (s ** 2).tolist()
    positive = s > 0
    if not np.any(positive):
        return np.zeros((xz.shape[0], 0)), eigvals, 0
    keep = positive & ((s[0] / np.maximum(s, 1e-12)) <= np.sqrt(cond_max))
    if not np.any(keep):
        keep[0] = True
    return u[:, keep], eigvals, int(np.count_nonzero(keep))


def save_nifti(data: np.ndarray, ref_img: nib.Nifti1Image, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    hdr = ref_img.header.copy()
    nib.save(nib.Nifti1Image(data, ref_img.affine, hdr), str(out_path))


def write_text_matrix(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, array, fmt="%.8f", delimiter="\t")


def run_acompcor(args: argparse.Namespace) -> Dict[str, object]:
    subdir = Path(args.subdir).resolve()
    in_path = Path(args.input).resolve()
    motion_path = Path(args.motion).resolve()
    censor_path = Path(args.censor).resolve() if args.censor else None
    out_base = Path(args.output_base).resolve()
    out_full_base = Path(args.output_full_base).resolve() if args.output_full_base else None

    in_img = load_nifti(in_path)
    brain_mask_path, ribbon_mask_path, aseg_path = choose_default_paths(subdir, args)
    brain_img = load_nifti(brain_mask_path)
    ribbon_img = load_nifti(ribbon_mask_path)
    aseg_img = load_nifti(aseg_path)

    data_2d, shape4d = load_data_matrix(in_img)
    n_vox, n_t = data_2d.shape

    brain_mask_3d = np.asanyarray(brain_img.dataobj) > 0
    ribbon_mask_3d = np.asanyarray(ribbon_img.dataobj) > 0
    if brain_mask_3d.shape != shape4d[:3]:
        raise RuntimeError(f"Brain mask shape mismatch: {brain_mask_3d.shape} vs input {shape4d[:3]}")
    if ribbon_mask_3d.shape != shape4d[:3]:
        raise RuntimeError(f"Ribbon mask shape mismatch: {ribbon_mask_3d.shape} vs input {shape4d[:3]}")

    motion = load_motion(motion_path, n_t)
    keep_mask = load_keep_mask(censor_path, n_t, args.censor_keep_value)
    if np.count_nonzero(keep_mask) < 2:
        raise RuntimeError("Too few uncensored frames to fit nuisance model.")

    wm_mask = clean_binary_mask(resample_label_mask(aseg_img, in_img, WM_LABELS), erode_iters=max(args.wm_erode_iters, 0))
    vent_mask = clean_binary_mask(resample_label_mask(aseg_img, in_img, VENT_LABELS))

    brain_dilated = clean_binary_mask(brain_mask_3d, dilate_iters=max(args.extraaxial_dilate_iters, 0))
    sd_pct = temporal_sd_percent(data_2d)
    extra_mask = (sd_pct.reshape(shape4d[:3]) > float(args.extraaxial_std_pct)) & (~brain_dilated)

    wm = compartment_regressors(data_2d, wm_mask, float(args.compartment_cond_max), "white_matter")
    vent = compartment_regressors(data_2d, vent_mask, float(args.compartment_cond_max), "ventricles")
    extra = compartment_regressors(data_2d, extra_mask, float(args.compartment_cond_max), "extra_axial")

    gs = compute_global_signal(data_2d, ribbon_mask_3d)
    gs_deriv = np.concatenate(([0.0], np.diff(gs)))

    design_pre = np.column_stack(
        [
            wm.regressors,
            vent.regressors,
            extra.regressors,
            gs[:, None],
            gs_deriv[:, None],
            motion,
        ]
    )
    design_basis, design_eigvals, n_design_components = condition_design(design_pre, float(args.design_cond_max))
    if design_basis.shape[1] == 0:
        raise RuntimeError("Combined nuisance design collapsed to zero columns.")

    brain_idx = np.flatnonzero(brain_mask_3d.reshape(-1) > 0)
    if brain_idx.size == 0:
        raise RuntimeError("Brain mask is empty.")

    y = data_2d[brain_idx, :].T
    intercept = np.ones((n_t, 1), dtype=np.float64)
    x_fit = np.column_stack([design_basis, intercept])
    betas, *_ = np.linalg.lstsq(x_fit[keep_mask, :], y[keep_mask, :], rcond=None)
    nuisance_only = design_basis @ betas[:-1, :]
    cleaned = y - nuisance_only
    cleaned_full = data_2d.copy()
    cleaned_full[brain_idx, :] = cleaned.T

    cleaned_full_4d = cleaned_full.reshape(shape4d).astype(np.float32)
    if keep_mask.all():
        cleaned_censored_4d = cleaned_full_4d
    else:
        cleaned_censored_4d = cleaned_full_4d[..., keep_mask].astype(np.float32)

    out_path = out_base.with_suffix("")
    if out_path.name.endswith(".nii"):
        out_path = out_path.with_suffix("")
    out_nii = out_path.with_suffix(".nii.gz")
    save_nifti(cleaned_censored_4d, in_img, out_nii)

    out_full_nii = None
    if out_full_base is not None:
        out_full_path = out_full_base.with_suffix("")
        if out_full_path.name.endswith(".nii"):
            out_full_path = out_full_path.with_suffix("")
        out_full_nii = out_full_path.with_suffix(".nii.gz")
        save_nifti(cleaned_full_4d, in_img, out_full_nii)

    mask_dir = out_nii.parent
    save_nifti(wm.mask.astype(np.uint8), brain_img, mask_dir / "WM_func_mask.nii.gz")
    save_nifti(vent.mask.astype(np.uint8), brain_img, mask_dir / "Ventricle_func_mask.nii.gz")
    save_nifti(extra.mask.astype(np.uint8), brain_img, mask_dir / "ExtraAxial_func_mask.nii.gz")

    write_text_matrix(mask_dir / "design_precondition.tsv", design_pre)
    write_text_matrix(mask_dir / "design_conditioned.tsv", design_basis)
    write_text_matrix(mask_dir / "motion.tsv", motion)
    write_text_matrix(mask_dir / "keep_mask.tsv", keep_mask.astype(np.int32)[:, None])

    summary = {
        "input": str(in_path),
        "output_censored": str(out_nii),
        "output_full": str(out_full_nii) if out_full_nii is not None else "",
        "brain_mask": str(brain_mask_path),
        "ribbon_mask": str(ribbon_mask_path),
        "aseg": str(aseg_path),
        "n_timepoints_in": int(n_t),
        "n_timepoints_kept": int(np.count_nonzero(keep_mask)),
        "n_timepoints_censored": int(n_t - np.count_nonzero(keep_mask)),
        "compartment_cond_max": float(args.compartment_cond_max),
        "design_cond_max": float(args.design_cond_max),
        "extraaxial_std_pct": float(args.extraaxial_std_pct),
        "extraaxial_dilate_iters": int(args.extraaxial_dilate_iters),
        "wm_erode_iters": int(args.wm_erode_iters),
        "wm": {
            "n_voxels": wm.n_voxels,
            "n_components": wm.n_components,
            "eigenvalues": wm.eigenvalues,
        },
        "ventricles": {
            "n_voxels": vent.n_voxels,
            "n_components": vent.n_components,
            "eigenvalues": vent.eigenvalues,
        },
        "extra_axial": {
            "n_voxels": extra.n_voxels,
            "n_components": extra.n_components,
            "eigenvalues": extra.eigenvalues,
        },
        "design_precondition_ncols": int(design_pre.shape[1]),
        "design_conditioned_ncols": int(n_design_components),
        "design_eigenvalues": design_eigvals,
    }
    (mask_dir / "acompcor_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    args = parse_args()
    summary = run_acompcor(args)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

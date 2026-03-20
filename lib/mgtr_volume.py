#!/usr/bin/env python3
"""MGTR volume denoising.

Replicates WorkingPipeline/res0urces/mgtr_volume.m behavior:
- Load in-brain mask (T1w_acpc_brain_func_mask).
- Load precomputed cortical ribbon mask in functional atlas space.
- Regress gray-matter global signal + intercept from each in-brain voxel.
- Add original voxel temporal mean back.
- Write:
  - Output_MGTR.nii.gz
  - Output_Betas.nii.gz
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import nibabel as nib
import numpy as np


def mgtr(subdir: str, input_nii: str, out_mgtr_base: str, out_betas_base: str) -> None:
    atlas_space = os.environ.get("AtlasSpace", "T1w")
    if atlas_space not in ("T1w", "MNINonlinear"):
        raise RuntimeError(f"Invalid AtlasSpace={atlas_space} (expected T1w or MNINonlinear)")

    xfms_dir = os.environ.get("FUNC_XFMS_DIRNAME", "rest")
    mask_path = os.path.join(subdir, "func", "xfms", xfms_dir, "T1w_acpc_brain_func_mask.nii.gz")
    if not os.path.isfile(mask_path):
        # Compatibility fallback for older pre-refactor outputs.
        mask_path = os.path.join(subdir, "func", "xfms", xfms_dir, "T1w_func_brain_mask.nii.gz")
    if not os.path.isfile(mask_path):
        raise FileNotFoundError(f"Missing mask: {mask_path}")

    mask_img = nib.load(mask_path)
    mask = np.asanyarray(mask_img.dataobj).reshape(-1)
    brain_idx = np.where(mask == 1)[0]
    if brain_idx.size == 0:
        raise RuntimeError(f"No in-brain voxels found in mask: {mask_path}")

    in_img = nib.load(input_nii)
    data = np.asanyarray(in_img.dataobj, dtype=np.float64)
    if data.ndim != 4:
        raise RuntimeError(f"Expected 4D input, got shape {data.shape}: {input_nii}")
    nx, ny, nz, nt = data.shape
    data2d = data.reshape((-1, nt))
    data_mean = data2d.mean(axis=1, keepdims=True)

    ribbon_source = os.environ.get("MGTR_RIBBON_SOURCE", "xfms").strip().lower()
    if ribbon_source not in ("xfms", "legacy_rois"):
        raise RuntimeError("MGTR_RIBBON_SOURCE must be 'xfms' or 'legacy_rois'")

    if ribbon_source == "legacy_rois":
        # Legacy MATLAB behavior: use func/rois/CorticalRibbon.nii.gz built from FreeSurfer lh/rh.ribbon.mgz.
        ribbon_path = os.path.join(subdir, "func", "rois", "CorticalRibbon.nii.gz")
        if not os.path.isfile(ribbon_path):
            subject = os.path.basename(os.path.normpath(subdir))
            t1_ref = mask_path
            lh_mgz = os.path.join(subdir, "anat", "T1w", subject, "mri", "lh.ribbon.mgz")
            rh_mgz = os.path.join(subdir, "anat", "T1w", subject, "mri", "rh.ribbon.mgz")
            lh_nii = os.path.join(subdir, "func", "rois", "lh.ribbon.nii.gz")
            rh_nii = os.path.join(subdir, "func", "rois", "rh.ribbon.nii.gz")
            os.makedirs(os.path.dirname(ribbon_path), exist_ok=True)
            if not (os.path.isfile(lh_mgz) and os.path.isfile(rh_mgz)):
                raise FileNotFoundError(
                    "Missing FreeSurfer ribbon mgz files required for MGTR_RIBBON_SOURCE=legacy_rois"
                )
            subprocess.run(
                ["mri_convert", "-i", lh_mgz, "-o", lh_nii, "--like", t1_ref],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["mri_convert", "-i", rh_mgz, "-o", rh_nii, "--like", t1_ref],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["fslmaths", lh_nii, "-add", rh_nii, ribbon_path],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["fslmaths", ribbon_path, "-bin", ribbon_path],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            for p in (lh_nii, rh_nii):
                try:
                    os.remove(p)
                except OSError:
                    pass
    else:
        ribbon_name = (
            "CorticalRibbon_nonlin_func_mask.nii.gz"
            if atlas_space == "MNINonlinear"
            else "CorticalRibbon_acpc_func_mask.nii.gz"
        )
        ribbon_path = os.path.join(subdir, "func", "xfms", xfms_dir, ribbon_name)
        if not os.path.isfile(ribbon_path):
            raise FileNotFoundError(
                f"Missing cortical ribbon mask in functional space for AtlasSpace={atlas_space}: {ribbon_path}"
            )
    gray = np.asanyarray(nib.load(ribbon_path).dataobj).reshape(-1)
    gray_idx = np.where(gray == 1)[0]
    if gray_idx.size == 0:
        raise RuntimeError(f"No ribbon voxels found in: {ribbon_path}")

    gs = data2d[gray_idx, :].mean(axis=0)
    x = np.column_stack((gs, np.ones(nt, dtype=np.float64)))
    x_pinv = np.linalg.pinv(x)

    y = data2d[brain_idx, :].T
    betas = x_pinv @ y
    resid = y - (x @ betas)
    data2d[brain_idx, :] = resid.T

    b = np.zeros((data2d.shape[0],), dtype=np.float64)
    b[brain_idx] = betas[0, :]

    out_data = (data2d + data_mean).reshape((nx, ny, nz, nt))
    out_b = b.reshape((nx, ny, nz))

    out_mgtr = out_mgtr_base + ".nii.gz"
    out_betas = out_betas_base + ".nii.gz"
    for p in (out_mgtr, out_betas):
        if os.path.exists(p):
            os.remove(p)

    nib.save(nib.Nifti1Image(out_data.astype(np.float64), in_img.affine, in_img.header), out_mgtr)
    nib.save(nib.Nifti1Image(out_b.astype(np.float64), mask_img.affine, mask_img.header), out_betas)


def main() -> int:
    ap = argparse.ArgumentParser(description="Run MGTR denoising")
    ap.add_argument("--subdir", required=True, help="Subject directory")
    ap.add_argument("--input", required=True, help="Input 4D timeseries (.nii.gz)")
    ap.add_argument("--output-mgtr-base", required=True, help="Output base path (without .nii.gz)")
    ap.add_argument("--output-betas-base", required=True, help="Output base path (without .nii.gz)")
    args = ap.parse_args()

    mgtr(args.subdir, args.input, args.output_mgtr_base, args.output_betas_base)
    return 0


if __name__ == "__main__":
    sys.exit(main())

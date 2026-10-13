#!/usr/bin/env python3
"""Spatial-feature ME-ICA reclassification.

This is a local/prototype successor to the NSI-only post-MEICA refinement.
It adds Moran's I and spectral slope through a model that predicts
log10(kappa/rho) from spatial information.

Design intent:

- NSI remains the canonical-network similarity term.
- A spatial model protects low-NSI but BOLD-like/noncanonical components.
- Kills require both low NSI and low predicted kappa/rho.
- Outputs are written to a new folder and do not modify the original tedana or
  Reclassify folders.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib
import numpy as np
import pandas as pd
import nibabel as nib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from nibabel.cifti2 import cifti2_axes
from nibabel.processing import resample_from_to
from scipy.io import loadmat
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nsi_tedana_model_pilot import (  # noqa: E402
    build_cortex_adjacency,
    extract_one_run,
    load_priors_fc,
    read_component_maps,
    ridge_r2_from_svd,
    spectral_basis,
    template_corr_details,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PRIORS = REPO_ROOT / "res0urces" / "Priors.mat"
DEFAULT_NEIGHBORS = REPO_ROOT / "lib" / "pfm-nsi" / "pfm_nsi" / "models" / "cifti_surf_neighbors_lr_normalwall.npz"
DEFAULT_TRAINING = REPO_ROOT / "res0urces" / "meica_reclassify_spatial_training_features.csv.gz"
WM_LABELS = np.array([2, 7, 41, 46, 77, 78, 79, 251, 252, 253, 254, 255], dtype=np.int32)
CEREBELLUM_LABELS = np.array([8, 47], dtype=np.int32)

SPATIAL_FEATURES = [
    "nsi",
    "morans_i",
    "spectral_slope",
]


def initial_keep_from_classification(df: pd.DataFrame) -> pd.Series:
    cls = df["classification"].astype(str).str.strip().str.lower()
    keep_labels = {"accepted", "provisionalaccept"}
    reject_labels = {"rejected", "ignored", "provisionalreject"}
    keep = cls.isin(keep_labels)
    unknown = ~(cls.isin(keep_labels) | cls.isin(reject_labels))
    if unknown.any():
        keep.loc[unknown] = False
    return keep


def train_spatial_model(training_csv: Path, args: argparse.Namespace):
    df = pd.read_csv(training_csv)
    missing = sorted(set(SPATIAL_FEATURES + ["log10_kappa_rho"]) - set(df.columns))
    if missing:
        raise ValueError(f"Training feature file missing required columns {missing}: {training_csv}")
    use = df[SPATIAL_FEATURES + ["log10_kappa_rho"]].copy()
    for col in use.columns:
        use[col] = pd.to_numeric(use[col], errors="coerce")
    use = use[np.isfinite(use["log10_kappa_rho"].to_numpy(dtype=float))].reset_index(drop=True)
    x = use[SPATIAL_FEATURES].to_numpy(dtype=float)
    y = use["log10_kappa_rho"].to_numpy(dtype=float)
    model = make_pipeline(
        SimpleImputer(strategy="median"),
        RandomForestRegressor(
            n_estimators=args.rf_trees,
            min_samples_leaf=args.rf_min_samples_leaf,
            n_jobs=-1,
            random_state=args.random_seed,
        ),
    )
    model.fit(x, y)
    return model


def compute_run_features(args: argparse.Namespace) -> pd.DataFrame:
    priors_fc = load_priors_fc(args.priors_mat)
    svd = np.linalg.svd(priors_fc, full_matrices=False)
    w = build_cortex_adjacency(priors_fc.shape[0], args.neighbor_table)
    spectral = spectral_basis(w, kmax=args.spectral_kmax)
    df = extract_one_run(
        ted_dir=args.tedana_dir,
        root=args.tedana_dir,
        priors_fc=priors_fc,
        svd=svd,
        ridge_lambda=args.ridge_lambda,
        w=w,
        spectral=spectral,
    )
    maps = read_component_maps(args.tedana_dir / "betas_OC.dtseries.nii", len(df), priors_fc.shape[0])
    maps_demeaned = maps - np.mean(maps, axis=0, keepdims=True)
    df["nsi_fc_demeaned"] = ridge_r2_from_svd(maps_demeaned, priors_fc, svd, lam=args.ridge_lambda)
    demeaned_corr, demeaned_signed_corr, demeaned_best_prior = template_corr_details(maps_demeaned, priors_fc)
    df["nsi_fc_demeaned_maxabs_corr"] = demeaned_corr
    df["nsi_fc_demeaned_signed_corr"] = demeaned_signed_corr
    df["nsi_fc_demeaned_best_prior"] = demeaned_best_prior
    df = add_tissue_loading_features(df, args)
    return df


def default_from_tedana(args: argparse.Namespace, name: str) -> Path:
    candidate = args.tedana_dir / name
    return candidate


def load_mask(path: Path, ref_img: nib.Nifti1Image, threshold: float = 0.5) -> np.ndarray:
    img = nib.load(str(path))
    if img.shape[:3] != ref_img.shape[:3] or not np.allclose(img.affine, ref_img.affine):
        img = resample_from_to(img, (ref_img.shape[:3], ref_img.affine), order=0)
    return np.asarray(img.dataobj) > threshold


def label_mask_from_aseg(path: Path, ref_img: nib.Nifti1Image, labels: np.ndarray) -> np.ndarray:
    img = nib.load(str(path))
    res = resample_from_to(img, (ref_img.shape[:3], ref_img.affine), order=0)
    data = np.asarray(res.dataobj).astype(np.int32)
    return np.isin(data, labels)


def mean_abs(data: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if not np.any(mask):
        return np.full(data.shape[3], np.nan, dtype=float)
    return np.mean(np.abs(data[mask, :]), axis=0)


def add_tissue_loading_features(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = df.copy()
    vol_path = args.component_nii if args.component_nii else args.tedana_dir / "betas_OC.nii.gz"
    if not vol_path.is_file():
        raise FileNotFoundError(f"Missing volumetric component map required for tissue features: {vol_path}")
    needed = [args.cortical_ribbon, args.subcortical_rois]
    if args.wm_mask is None and args.aseg is None:
        raise FileNotFoundError("Tissue features require either --wm-mask or --aseg")
    if any(p is None or not Path(p).is_file() for p in needed):
        missing = [str(p) for p in needed if p is None or not Path(p).is_file()]
        raise FileNotFoundError(f"Missing tissue mask(s) required for reclassification: {missing}")

    img = nib.load(str(vol_path))
    data = np.asarray(img.dataobj, dtype=np.float32)
    if data.ndim != 4 or data.shape[3] != len(out):
        raise ValueError(
            f"Expected 4D component NIfTI with {len(out)} components, got shape {data.shape}: {vol_path}"
        )

    if args.wm_mask is not None and args.wm_mask.is_file():
        wm_mask = load_mask(args.wm_mask, img, threshold=0.5)
    else:
        wm_mask = label_mask_from_aseg(args.aseg, img, WM_LABELS)
    cortex_mask = load_mask(args.cortical_ribbon, img, threshold=0.5)
    subcort_mask = load_mask(args.subcortical_rois, img, threshold=0.5)
    cerebellum_mask = label_mask_from_aseg(args.subcortical_rois, img, CEREBELLUM_LABELS)
    gray_mask = cortex_mask | subcort_mask
    required_masks = {
        "white_matter": wm_mask,
        "cortical_ribbon": cortex_mask,
        "subcortical_rois": subcort_mask,
        "cerebellum": cerebellum_mask,
    }
    empty = [name for name, mask in required_masks.items() if not np.any(mask)]
    if empty:
        raise ValueError(f"Required tissue mask(s) have zero voxels after resampling: {empty}")

    wm = mean_abs(data, wm_mask)
    cortex = mean_abs(data, cortex_mask)
    subcort = mean_abs(data, subcort_mask)
    cerebellum = mean_abs(data, cerebellum_mask)
    gray = mean_abs(data, gray_mask)
    eps = np.finfo(float).eps
    out["tissue_metrics_available"] = True
    out["wm_abs_mean"] = wm
    out["cortex_abs_mean"] = cortex
    out["subcort_abs_mean"] = subcort
    out["cerebellum_abs_mean"] = cerebellum
    out["gray_abs_mean"] = gray
    out["wm_to_gray_abs_mean"] = wm / np.maximum(gray, eps)
    out["wm_to_cortex_abs_mean"] = wm / np.maximum(cortex, eps)
    out["subcort_to_cortex_abs_mean"] = subcort / np.maximum(cortex, eps)
    out["cerebellum_to_cortex_abs_mean"] = cerebellum / np.maximum(cortex, eps)
    out["wm_mask_voxels"] = int(wm_mask.sum())
    out["cortex_mask_voxels"] = int(cortex_mask.sum())
    out["subcort_mask_voxels"] = int(subcort_mask.sum())
    out["cerebellum_mask_voxels"] = int(cerebellum_mask.sum())
    return out


def adaptive_nsi_threshold(logratio: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    if args.kill_mode == "adaptive":
        return np.clip(
            args.kill_intercept - args.kill_slope * logratio,
            args.kill_nsi_min,
            args.kill_nsi_max,
        )
    return np.full_like(logratio, args.kill_nsi, dtype=float)


def refine(df: pd.DataFrame, model, args: argparse.Namespace) -> pd.DataFrame:
    out = df.copy()
    numeric_cols = ["kappa", "rho", "nsi", "nsi_fc_demeaned", "log10_kappa_rho"] + SPATIAL_FEATURES
    tissue_cols = [
        "wm_to_gray_abs_mean",
        "wm_to_cortex_abs_mean",
        "subcort_to_cortex_abs_mean",
        "cerebellum_to_cortex_abs_mean",
    ]
    for col in numeric_cols + [c for c in tissue_cols if c in out.columns]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    x = out[SPATIAL_FEATURES].to_numpy(dtype=float)
    out["pred_log10_kappa_rho_spatial"] = model.predict(x)
    out["pred_residual_log10_kappa_rho"] = out["log10_kappa_rho"] - out["pred_log10_kappa_rho_spatial"]

    keep_init = initial_keep_from_classification(out)
    if args.guardrail_kappa_rho:
        keep_init &= out["kappa"] > out["rho"]
    out["keep_init"] = keep_init

    logratio = out["log10_kappa_rho"].to_numpy(dtype=float)
    pred = out["pred_log10_kappa_rho_spatial"].to_numpy(dtype=float)
    decision_nsi = out[args.decision_nsi_col].to_numpy(dtype=float)
    raw_nsi = out["nsi"].to_numpy(dtype=float)
    nsi_thresh = adaptive_nsi_threshold(logratio, args)
    out["adaptive_nsi_kill_threshold"] = nsi_thresh
    out["decision_nsi"] = decision_nsi

    # Candidate classes. Low NSI but high predicted logratio is not auto-killed:
    # this is the BOLD-like/noncanonical review zone.
    low_nsi = decision_nsi < nsi_thresh
    low_pred = pred < args.pred_kill_logratio
    high_pred = pred >= args.pred_protect_logratio
    high_nsi = decision_nsi >= args.rescue_nsi
    observed_bold_like = logratio >= args.observed_protect_logratio
    high_morans = out["morans_i"].to_numpy(dtype=float) >= args.protect_morans
    smooth_slope = out["spectral_slope"].to_numpy(dtype=float) <= args.protect_spectral_slope_max
    high_demeaned_nsi = out["nsi_fc_demeaned"].to_numpy(dtype=float) >= args.protect_demeaned_nsi
    high_template = decision_nsi >= args.protect_template_corr
    old_spatial_protect = high_pred | observed_bold_like | high_template | high_morans | smooth_slope | high_demeaned_nsi

    tissue_available = out.get("tissue_metrics_available", pd.Series(False, index=out.index)).astype(bool).to_numpy()
    if "wm_to_gray_abs_mean" in out.columns:
        wm_to_gray = out["wm_to_gray_abs_mean"].to_numpy(dtype=float)
        sub_to_cortex = out["subcort_to_cortex_abs_mean"].to_numpy(dtype=float)
        cereb_to_cortex = out["cerebellum_to_cortex_abs_mean"].to_numpy(dtype=float)
    else:
        wm_to_gray = np.full(len(out), np.nan)
        sub_to_cortex = np.full(len(out), np.nan)
        cereb_to_cortex = np.full(len(out), np.nan)

    # Use WM relative to cortical ribbon. WM:gray can be diluted when a vascular
    # component also loads strongly in subcortical/cerebellar tissue.
    likely_wm_vascular = tissue_available & (
        out["wm_to_cortex_abs_mean"].to_numpy(dtype=float) >= args.wm_kill_ratio
    )
    subcortical_or_cerebellar = tissue_available & (
        (sub_to_cortex >= args.subcort_protect_ratio) | (cereb_to_cortex >= args.cerebellum_protect_ratio)
    )
    protect_subcortical_review = (
        keep_init.to_numpy(dtype=bool)
        & low_nsi
        & subcortical_or_cerebellar
        & ~likely_wm_vascular
    )
    review_low_nsi_bold = keep_init.to_numpy(dtype=bool) & low_nsi & old_spatial_protect & ~protect_subcortical_review

    kill_candidate = keep_init.to_numpy(dtype=bool) & low_nsi & ~protect_subcortical_review
    protected_low_nsi_bold = protect_subcortical_review

    rescue_candidate = (~keep_init.to_numpy(dtype=bool)) & high_nsi & (pred >= args.pred_rescue_logratio)
    if args.guardrail_kappa_rho:
        rescue_candidate &= out["kappa"].to_numpy(dtype=float) > out["rho"].to_numpy(dtype=float)

    # Optional cap/prioritization: kill the strongest low-NSI + low-pred evidence.
    killed = kill_candidate.copy()
    if args.kill_priority_enable and np.any(kill_candidate):
        idx = np.where(kill_candidate)[0]
        vcol = "variance explained" if "variance explained" in out.columns else "normalized variance explained"
        if vcol in out.columns:
            var = pd.to_numeric(out[vcol], errors="coerce").fillna(0).to_numpy(dtype=float)
        else:
            var = np.ones(len(out), dtype=float)

        def z(a):
            a = np.asarray(a, dtype=float)
            sd = np.nanstd(a)
            if sd <= 0:
                return np.zeros_like(a)
            return (a - np.nanmean(a)) / sd

        evidence = (
            args.kill_priority_w_nsi * z(-(nsi_thresh[idx] - decision_nsi[idx]))
            + args.kill_priority_w_pred * z(-(pred[idx]))
            + args.kill_priority_w_var * z(var[idx])
        )
        order = idx[np.argsort(-evidence)]
        if args.kill_max_count > 0:
            order = order[: args.kill_max_count]
        killed[:] = False
        killed[order] = True

    keep_final = keep_init.to_numpy(dtype=bool)
    keep_final[rescue_candidate] = True
    keep_final[killed] = False
    if args.guardrail_kappa_rho:
        keep_final &= out["kappa"].to_numpy(dtype=float) > out["rho"].to_numpy(dtype=float)

    out["rescue_candidate_spatial"] = rescue_candidate
    out["kill_candidate_spatial"] = kill_candidate
    out["protect_high_pred"] = high_pred
    out["protect_observed_bold_like"] = observed_bold_like
    out["protect_high_template_corr"] = high_template
    out["protect_high_morans_i"] = high_morans
    out["protect_smooth_spectral_slope"] = smooth_slope
    out["protect_high_demeaned_nsi"] = high_demeaned_nsi
    out["would_old_spatial_protect"] = old_spatial_protect
    out["likely_wm_vascular"] = likely_wm_vascular
    out["subcortical_or_cerebellar_loading"] = subcortical_or_cerebellar
    out["protected_subcortical_or_cerebellar_review"] = protect_subcortical_review
    out["review_low_nsi_bold_like"] = review_low_nsi_bold
    out["raw_nsi_clipped_zero"] = raw_nsi == 0
    out["protected_low_nsi_bold_like"] = protected_low_nsi_bold
    out["rescued"] = rescue_candidate
    out["killed"] = killed
    out["accepted_final"] = keep_final

    category = np.full(len(out), "neutral", dtype=object)
    category[keep_init.to_numpy(dtype=bool)] = "init_accept"
    category[~keep_init.to_numpy(dtype=bool)] = "init_reject"
    category[review_low_nsi_bold] = "killed_low_nsi_bold_like_review"
    category[likely_wm_vascular & keep_init.to_numpy(dtype=bool) & low_nsi] = "killed_likely_wm_vascular"
    category[protected_low_nsi_bold] = "protected_subcortical_or_cerebellar_review"
    category[rescue_candidate] = "rescued"
    category[killed] = "killed"
    category[killed & likely_wm_vascular] = "killed_likely_wm_vascular"
    category[killed & review_low_nsi_bold & ~likely_wm_vascular] = "killed_low_nsi_bold_like_review"
    category[protected_low_nsi_bold] = "protected_subcortical_or_cerebellar_review"
    out["spatial_decision_category"] = category
    return out


def write_ids(path: Path, ids: List[int]) -> None:
    path.write_text("\n".join(str(int(i)) for i in ids) + ("\n" if ids else ""), encoding="utf-8")


def make_traffic_plot(out: pd.DataFrame, out_png: Path, args: argparse.Namespace) -> None:
    x = out["log10_kappa_rho"].to_numpy(dtype=float)
    y = out["decision_nsi"].to_numpy(dtype=float)
    pred = out["pred_log10_kappa_rho_spatial"].to_numpy(dtype=float)
    keep_init = out["keep_init"].to_numpy(dtype=bool)
    killed = out["killed"].to_numpy(dtype=bool)
    rescued = out["rescued"].to_numpy(dtype=bool)
    protected = out["protected_low_nsi_bold_like"].to_numpy(dtype=bool)

    xmin = min(np.nanmin(x) - 0.05, -0.5)
    xmax = max(np.nanmax(x) + 0.05, 1.2)
    ymax = max(np.nanmax(y) + 0.05, 0.65)
    xline = np.linspace(xmin, xmax, 400)
    thresh = adaptive_nsi_threshold(xline, args)

    plt.figure(figsize=(8.2, 5.2), dpi=160)
    ax = plt.gca()
    ax.axvspan(xmin, 0.0, color="#f4cccc", alpha=0.35, label="low observed kappa/rho")
    ax.fill_between(xline, 0, thresh, color="#f4cccc", alpha=0.55, label="kill zone: low NSI")
    ax.fill_between(xline, thresh, ymax, where=xline >= args.observed_protect_logratio, color="#d9ead3", alpha=0.35, label="protect zone")
    ax.fill_between(xline, thresh, ymax, where=(xline >= 0) & (xline < args.observed_protect_logratio), color="#fff2cc", alpha=0.35, label="review zone")
    ax.plot(xline, thresh, "k--", linewidth=1.2, label="adaptive NSI threshold")
    finite_pred = pred[np.isfinite(pred)]
    vmin = float(np.nanmin(finite_pred)) if finite_pred.size else None
    vmax = float(np.nanmax(finite_pred)) if finite_pred.size else None

    scatter = ax.scatter(
        x[~keep_init],
        y[~keep_init],
        c=pred[~keep_init],
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        s=20,
        alpha=0.65,
        marker="o",
        edgecolors="none",
        label="init reject",
    )
    ax.scatter(
        x[keep_init],
        y[keep_init],
        c=pred[keep_init],
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        s=24,
        alpha=0.75,
        marker="o",
        edgecolors="#444444",
        linewidths=0.25,
        label="init accept",
    )
    if np.any(protected):
        ax.scatter(
            x[protected],
            y[protected],
            s=72,
            facecolors="none",
            edgecolors="#2ca02c",
            linewidths=1.2,
            label="Marked for review",
        )
    if np.any(rescued):
        ax.scatter(x[rescued], y[rescued], s=76, c="#ffd92f", edgecolors="black", linewidths=0.5, label="rescued")
    if np.any(killed):
        ax.scatter(x[killed], y[killed], s=76, c="#00d5d8", edgecolors="black", linewidths=0.5, label="killed")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.02, ymax)
    ax.set_xlabel("observed log10(kappa / rho)")
    ax.set_ylabel("NSI")
    ax.set_title("ME-ICA reclassification based on spatial properties")
    ax.grid(color="#dddddd", linewidth=0.5)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("spatial model predicted log10(kappa/rho)")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def make_predicted_observed_plot(out: pd.DataFrame, out_png: Path) -> None:
    x = out["log10_kappa_rho"].to_numpy(dtype=float)
    y = out["pred_log10_kappa_rho_spatial"].to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    r = np.corrcoef(x[keep], y[keep])[0, 1] if keep.sum() > 2 else np.nan
    plt.figure(figsize=(5.0, 4.5), dpi=160)
    color = out["decision_nsi"] if "decision_nsi" in out.columns else out["nsi"]
    plt.scatter(x, y, s=24, c=color, cmap="viridis", alpha=0.8, edgecolors="none")
    lo = min(np.nanmin(x), np.nanmin(y))
    hi = max(np.nanmax(x), np.nanmax(y))
    plt.plot([lo, hi], [lo, hi], color="black", linewidth=1)
    plt.xlabel("observed log10(kappa/rho)")
    plt.ylabel("predicted from spatial features")
    plt.title(f"Spatial model agreement on this run (r={r:.2f})")
    cb = plt.colorbar()
    cb.set_label("decision NSI")
    plt.grid(color="#dddddd", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def summarize(out: pd.DataFrame, args: argparse.Namespace) -> Dict[str, object]:
    init_accept = int(out["keep_init"].sum())
    final_accept = int(out["accepted_final"].sum())
    killed = int(out["killed"].sum())
    rescued = int(out["rescued"].sum())
    protected = int(out["protected_low_nsi_bold_like"].sum())
    likely_wm = int(out.get("likely_wm_vascular", pd.Series(False, index=out.index)).sum())
    subcort_review = int(out.get("protected_subcortical_or_cerebellar_review", pd.Series(False, index=out.index)).sum())
    low_nsi_review = int(out.get("review_low_nsi_bold_like", pd.Series(False, index=out.index)).sum())
    return {
        "tedana_dir": str(args.tedana_dir),
        "training_features": str(args.training_features),
        "n_components": int(len(out)),
        "init_accepted": init_accept,
        "final_accepted": final_accept,
        "killed": killed,
        "rescued": rescued,
        "protected_low_nsi_bold_like": protected,
        "protected_subcortical_or_cerebellar_review": subcort_review,
        "likely_wm_vascular": likely_wm,
        "review_low_nsi_bold_like": low_nsi_review,
        "settings": {
            "decision_nsi_col": args.decision_nsi_col,
            "kill_mode": args.kill_mode,
            "kill_nsi": args.kill_nsi,
            "kill_nsi_min": args.kill_nsi_min,
            "kill_nsi_max": args.kill_nsi_max,
            "kill_intercept": args.kill_intercept,
            "kill_slope": args.kill_slope,
            "pred_kill_logratio": args.pred_kill_logratio,
            "pred_protect_logratio": args.pred_protect_logratio,
            "observed_protect_logratio": args.observed_protect_logratio,
            "rescue_nsi": args.rescue_nsi,
            "pred_rescue_logratio": args.pred_rescue_logratio,
            "protect_template_corr": args.protect_template_corr,
            "protect_morans": args.protect_morans,
            "protect_spectral_slope_max": args.protect_spectral_slope_max,
            "protect_demeaned_nsi": args.protect_demeaned_nsi,
            "wm_kill_ratio": args.wm_kill_ratio,
            "subcort_protect_ratio": args.subcort_protect_ratio,
            "cerebellum_protect_ratio": args.cerebellum_protect_ratio,
        },
    }


def class_masks(out: pd.DataFrame) -> Dict[str, pd.Series]:
    return {
        "AcceptedFinal": out["accepted_final"].astype(bool),
        "RejectedFinal": ~out["accepted_final"].astype(bool),
        "Killed": out["killed"].astype(bool),
        "Rescued": out["rescued"].astype(bool),
        "AcceptedMarkedForReview": out["protected_subcortical_or_cerebellar_review"].astype(bool),
        "LowNSIBoldLikeReview": out["review_low_nsi_bold_like"].astype(bool),
        "LikelyWMVascular": out["likely_wm_vascular"].astype(bool),
    }


def write_qc_ids(path: Path, ids: List[int]) -> None:
    write_ids(path, ids)


def write_nifti_subset(img: nib.Nifti1Image, ids: List[int], out_path: Path) -> None:
    data = np.asarray(img.dataobj, dtype=np.float32)
    if data.ndim != 4:
        raise ValueError(f"Expected 4D component NIfTI, got shape {data.shape}")
    valid = [i for i in ids if 0 <= i < data.shape[3]]
    if valid:
        subset = data[:, :, :, valid]
    else:
        subset = np.zeros(data.shape[:3] + (1,), dtype=np.float32)
    hdr = img.header.copy()
    hdr.set_data_shape(subset.shape)
    nib.save(nib.Nifti1Image(subset, img.affine, hdr), str(out_path))


def write_cifti_subset(cifti_path: Path, ids: List[int], out_path: Path) -> None:
    if not cifti_path.is_file():
        return
    img = nib.load(str(cifti_path))
    data = np.asarray(img.dataobj, dtype=np.float32)
    ax0 = img.header.get_axis(0)
    ax1 = img.header.get_axis(1)
    if isinstance(ax1, cifti2_axes.BrainModelAxis):
        brain_axis = ax1
        valid = [i for i in ids if 0 <= i < data.shape[0]]
        maps = data[valid, :] if valid else np.zeros((1, data.shape[1]), dtype=np.float32)
    elif isinstance(ax0, cifti2_axes.BrainModelAxis):
        brain_axis = ax0
        valid = [i for i in ids if 0 <= i < data.shape[1]]
        maps = data[:, valid].T if valid else np.zeros((1, data.shape[0]), dtype=np.float32)
    else:
        raise ValueError(f"Could not find BrainModelAxis in CIFTI: {cifti_path}")
    names = [f"ica_{i:03d}" for i in valid] if valid else ["empty"]
    scalar_axis = cifti2_axes.ScalarAxis(names)
    header = nib.cifti2.Cifti2Header.from_axes((scalar_axis, brain_axis))
    nib.save(nib.Cifti2Image(maps, header=header, nifti_header=img.nifti_header), str(out_path))


def resample_bg_to_component(bg_path: Path | None, ref_img: nib.Nifti1Image) -> np.ndarray | None:
    if bg_path is None or not bg_path.is_file():
        return None
    bg_img = nib.load(str(bg_path))
    if bg_img.shape[:3] != ref_img.shape[:3] or not np.allclose(bg_img.affine, ref_img.affine):
        bg_img = resample_from_to(bg_img, (ref_img.shape[:3], ref_img.affine), order=1)
    bg = np.asarray(bg_img.dataobj, dtype=np.float32)
    if bg.ndim == 4:
        bg = bg[:, :, :, 0]
    finite = np.isfinite(bg)
    if not np.any(finite):
        return None
    lo, hi = np.percentile(bg[finite], [1, 99])
    if hi <= lo:
        return np.zeros_like(bg, dtype=np.float32)
    return np.clip((bg - lo) / (hi - lo), 0, 1)


def slice_indices_for_component(comp: np.ndarray, n_slices: int = 8) -> Dict[int, np.ndarray]:
    abs_comp = np.abs(comp)
    finite = np.isfinite(abs_comp)
    if not np.any(finite):
        return {axis: np.linspace(0, comp.shape[axis] - 1, n_slices).astype(int) for axis in range(3)}
    thresh = np.nanpercentile(abs_comp[finite], 97.5)
    coords = np.array(np.where(abs_comp >= max(thresh, np.finfo(float).eps)))
    if coords.size == 0:
        coords = np.array(np.where(finite))
    out: Dict[int, np.ndarray] = {}
    for axis in range(3):
        vals = coords[axis]
        if vals.size == 0:
            lo, hi = 0, comp.shape[axis] - 1
        else:
            lo, hi = np.percentile(vals, [5, 95])
            pad = max(2, 0.05 * comp.shape[axis])
            lo = max(0, int(math.floor(lo - pad)))
            hi = min(comp.shape[axis] - 1, int(math.ceil(hi + pad)))
            if hi <= lo:
                center = int(np.median(vals))
                half = max(3, n_slices // 2)
                lo = max(0, center - half)
                hi = min(comp.shape[axis] - 1, center + half)
        idx = np.unique(np.linspace(lo, hi, n_slices).round().astype(int))
        if idx.size < n_slices:
            idx = np.unique(np.linspace(0, comp.shape[axis] - 1, n_slices).round().astype(int))
        out[axis] = idx[:n_slices]
    return out


def load_component_timecourses(tedana_dir: Path, n_components: int) -> np.ndarray | None:
    for name in ("desc-ICA_mixing.tsv", "ica_mixing.tsv"):
        path = tedana_dir / name
        if not path.is_file():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] < n_components:
            continue
        return df.iloc[:, :n_components].to_numpy(dtype=float)
    return None


def load_motion_params(path: Path | None) -> np.ndarray | None:
    if path is None or not path.is_file():
        return None
    motion = np.loadtxt(path, dtype=float)
    if motion.ndim == 1:
        motion = motion.reshape(1, -1)
    if motion.shape[1] < 6:
        return None
    return motion[:, :6]


def compute_fd_from_motion(motion: np.ndarray) -> np.ndarray:
    diff = np.zeros_like(motion, dtype=float)
    diff[1:, :] = np.abs(np.diff(motion, axis=0))
    # FSL MCFLIRT convention is rotations first, then translations. Convert
    # rotational displacement to mm on a 50 mm radius sphere.
    return np.sum(diff[:, 3:6], axis=1) + 50.0 * np.sum(diff[:, 0:3], axis=1)


def corr_1d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = min(a.size, b.size)
    if n < 3:
        return np.nan
    a = a[:n] - np.nanmean(a[:n])
    b = b[:n] - np.nanmean(b[:n])
    den = float(np.sqrt(np.nansum(a * a) * np.nansum(b * b)))
    if den <= 0:
        return np.nan
    return float(np.nansum(a * b) / den)


def add_motion_coupling_features(out: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    tc = load_component_timecourses(args.tedana_dir, len(out))
    motion = load_motion_params(args.motion_params)
    if tc is None or motion is None:
        return out
    fd = compute_fd_from_motion(motion)
    labels = ["rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z"]
    max_abs = np.full(len(out), np.nan, dtype=float)
    fd_corr = np.full(len(out), np.nan, dtype=float)
    best_label: List[str] = [""] * len(out)
    for pos, row in out.iterrows():
        cid = int(row["component_id"])
        if cid < 0 or cid >= tc.shape[1]:
            continue
        y = tc[:, cid]
        corrs = np.array([corr_1d(y, motion[:, j]) for j in range(6)], dtype=float)
        abs_corrs = np.abs(corrs)
        if np.any(np.isfinite(abs_corrs)):
            best = int(np.nanargmax(abs_corrs))
            max_abs[pos] = float(abs_corrs[best])
            best_label[pos] = labels[best]
        fd_corr[pos] = abs(corr_1d(y, fd))
    out = out.copy()
    out["motion_max_abs_corr"] = max_abs
    out["motion_best_param"] = best_label
    out["motion_fd_abs_corr"] = fd_corr
    out["motion_artifact_like"] = (max_abs >= args.motion_corr_threshold) | (fd_corr >= args.motion_corr_threshold)
    return out


def add_timecourse_spectral_features(out: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    tc = load_component_timecourses(args.tedana_dir, len(out))
    if tc is None:
        return out
    tr = args.tr
    if tr is None and args.tr_file is not None and args.tr_file.is_file():
        vals = args.tr_file.read_text().split()
        if vals:
            tr = float(vals[0])
    if tr is None or tr <= 0:
        tr = 1.0
    n_t = tc.shape[0]
    freq = np.fft.rfftfreq(n_t, d=float(tr))
    window = np.hanning(n_t)
    peak_freq = np.full(len(out), np.nan, dtype=float)
    peak_power = np.full(len(out), np.nan, dtype=float)
    band_003_007 = np.full(len(out), np.nan, dtype=float)
    band_0045_006 = np.full(len(out), np.nan, dtype=float)
    total_low = np.full(len(out), np.nan, dtype=float)
    for pos, row in out.iterrows():
        cid = int(row["component_id"])
        if cid < 0 or cid >= tc.shape[1]:
            continue
        y = np.asarray(tc[:, cid], dtype=float)
        y = y - np.nanmean(y)
        if n_t > 2:
            x = np.arange(n_t, dtype=float)
            coef = np.polyfit(x, y, 1)
            y = y - np.polyval(coef, x)
        power = np.abs(np.fft.rfft(y * window)) ** 2
        low = (freq > 0) & (freq <= args.tc_peak_max_hz)
        if np.any(low):
            low_idx = np.where(low)[0]
            best = low_idx[int(np.argmax(power[low_idx]))]
            peak_freq[pos] = float(freq[best])
            peak_power[pos] = float(power[best])
            total_low[pos] = float(np.sum(power[low_idx]))
        m = (freq >= 0.03) & (freq <= 0.07)
        if np.any(m):
            band_003_007[pos] = float(np.sum(power[m]))
        m = (freq >= 0.045) & (freq <= 0.060)
        if np.any(m):
            band_0045_006[pos] = float(np.sum(power[m]))
    out = out.copy()
    out["tc_tr_seconds"] = float(tr)
    out["tc_peak_freq_hz_le_0p2"] = peak_freq
    out["tc_peak_power_le_0p2"] = peak_power
    out["tc_power_0p03_0p07_hz"] = band_003_007
    out["tc_power_0p045_0p060_hz"] = band_0045_006
    out["tc_power_le_0p2_hz"] = total_low
    out["tc_frac_power_0p03_0p07_hz"] = band_003_007 / np.maximum(total_low, np.finfo(float).eps)
    out["tc_frac_power_0p045_0p060_hz"] = band_0045_006 / np.maximum(total_low, np.finfo(float).eps)
    return out


def load_prior_network_labels(priors_mat: Path) -> List[str]:
    try:
        mat = loadmat(str(priors_mat))
        labels = mat["Priors"]["NetworkLabels"][0, 0]
    except Exception:
        return []
    out: List[str] = []
    for item in np.ravel(labels):
        try:
            value = str(np.ravel(item)[0])
        except Exception:
            value = str(item)
        out.append(value)
    return out


def network_family(label: str) -> str:
    if not label:
        return "Network"
    family = label.split("_", 1)[0].split("/", 1)[0]
    if family == "SomatoCognitiveAction":
        return "SCA"
    return family


def component_title_class(row: pd.Series, prior_labels: List[str]) -> str:
    if bool(row.get("accepted_final", False)):
        prior = row.get("nsi_fc_demeaned_best_prior", row.get("nsi_fc_best_prior", np.nan))
        label = ""
        if np.isfinite(prior):
            idx = int(prior) - 1
            if 0 <= idx < len(prior_labels):
                label = network_family(prior_labels[idx])
        return f"Accepted: {label}-like" if label else "Accepted"
    if bool(row.get("likely_wm_vascular", False)):
        return "Artifact: WM/vascular"
    if bool(row.get("motion_artifact_like", False)):
        return "Artifact: motion"
    return "Artifact"


def component_spectrum(
    timecourse: np.ndarray | None,
    tr: float | None,
) -> tuple[str, np.ndarray | None, np.ndarray | None, float | None]:
    if timecourse is None or timecourse.size == 0:
        return "peak unavailable", None, None, None
    tr = float(tr) if tr and tr > 0 else 1.0
    y = np.asarray(timecourse, dtype=float)
    y = y - np.nanmean(y)
    if y.size > 2:
        x = np.arange(y.size, dtype=float)
        y = y - np.polyval(np.polyfit(x, y, 1), x)
    freq = np.fft.rfftfreq(y.size, d=tr)
    power = np.abs(np.fft.rfft(y * np.hanning(y.size))) ** 2
    mask = (freq > 0) & (freq <= 0.20)
    if not np.any(mask):
        return "peak unavailable", freq, power, None
    idx = np.where(mask)[0]
    peak_idx = idx[int(np.argmax(power[idx]))]
    peak = float(freq[peak_idx])
    return f"peak {peak:.3f} Hz", freq, power, peak


def numeric_value(row: pd.Series, name: str) -> float:
    try:
        value = float(row.get(name, np.nan))
    except Exception:
        value = np.nan
    return value if np.isfinite(value) else np.nan


def fmt_value(value: float, digits: int = 3, suffix: str = "") -> str:
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}{suffix}"


def variance_scale_label(variance_explained: float) -> str:
    if not np.isfinite(variance_explained):
        return "unknown"
    if variance_explained >= 1.0:
        return "large"
    if variance_explained >= 0.1:
        return "small"
    return "tiny"


def setup_readout_axis(ax, title: str, accent: str = "#333333") -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#cccccc")
        spine.set_linewidth(0.6)
    ax.set_facecolor("#fbfbf8")
    ax.text(0.04, 0.90, title, ha="left", va="top", fontsize=7.5, fontweight="bold", color=accent)


def add_readout_bar(
    ax,
    y: float,
    label: str,
    value: float,
    vmax: float,
    color: str,
    text: str,
    threshold: float | None = None,
) -> None:
    ax.text(0.04, y + 0.035, label, ha="left", va="bottom", fontsize=6.7, color="#333333")
    ax.text(0.96, y + 0.035, text, ha="right", va="bottom", fontsize=6.7, color="#333333")
    ax.add_patch(plt.Rectangle((0.04, y - 0.005), 0.92, 0.035, color="#e8e8e4", linewidth=0))
    if np.isfinite(value):
        frac = float(np.clip(value / max(vmax, np.finfo(float).eps), 0.0, 1.0))
        ax.add_patch(plt.Rectangle((0.04, y - 0.005), 0.92 * frac, 0.035, color=color, linewidth=0))
    if threshold is not None and np.isfinite(threshold):
        tx = 0.04 + 0.92 * float(np.clip(threshold / max(vmax, np.finfo(float).eps), 0.0, 1.0))
        ax.plot([tx, tx], [y - 0.015, y + 0.04], color="#222222", linewidth=0.7, linestyle="--")


def add_decision_readout(ax, row: pd.Series) -> None:
    setup_readout_axis(ax, "Decision / variance", "#2c5aa0")
    nsi = numeric_value(row, "decision_nsi")
    threshold = numeric_value(row, "adaptive_nsi_kill_threshold")
    ve = numeric_value(row, "variance explained")
    nve = numeric_value(row, "normalized variance explained")
    add_readout_bar(
        ax,
        0.68,
        "NSI",
        nsi,
        0.30,
        "#4c78a8",
        f"{fmt_value(nsi)}  thr {fmt_value(threshold)}",
        threshold=threshold,
    )
    add_readout_bar(
        ax,
        0.47,
        "Variance explained",
        ve,
        2.0,
        "#72b7b2",
        f"{fmt_value(ve, 3, '%')}  {variance_scale_label(ve)}",
    )
    add_readout_bar(
        ax,
        0.26,
        "Norm variance",
        nve,
        2.0,
        "#54a24b",
        fmt_value(nve, 3, "%"),
    )
    category = str(row.get("spatial_decision_category", ""))
    ax.text(0.04, 0.09, category.replace("_", " "), ha="left", va="bottom", fontsize=6.6, color="#555555")


def add_model_fit_readout(ax, row: pd.Series, model: str) -> None:
    if model == "rho":
        setup_readout_axis(ax, "Rho / S0", "#8c6d31")
        score = numeric_value(row, "rho")
        dice = numeric_value(row, "dice_FS0")
        count = numeric_value(row, "countsigFS0")
        color = "#b279a2"
        fit_label = "S0 fit"
    else:
        setup_readout_axis(ax, "Kappa / T2*", "#1b6a68")
        score = numeric_value(row, "kappa")
        dice = numeric_value(row, "dice_FT2")
        count = numeric_value(row, "countsigFT2")
        color = "#59a14f"
        fit_label = "T2* fit"
    add_readout_bar(ax, 0.60, model.capitalize(), score, 60.0, color, fmt_value(score, 2))
    add_readout_bar(ax, 0.36, fit_label, dice, 0.75, color, fmt_value(dice, 3))
    count_text = f"{int(count):d} vox" if np.isfinite(count) else "n/a"
    ax.text(0.04, 0.13, count_text, ha="left", va="bottom", fontsize=6.7, color="#333333")


def add_spatial_context_readout(ax, row: pd.Series, peak_text: str) -> None:
    setup_readout_axis(ax, "Spatial / artifact context", "#6f4e7c")
    wm_ctx = numeric_value(row, "wm_to_cortex_abs_mean")
    sub_ctx = numeric_value(row, "subcort_to_cortex_abs_mean")
    cereb_ctx = numeric_value(row, "cerebellum_to_cortex_abs_mean")
    motion = numeric_value(row, "motion_max_abs_corr")
    add_readout_bar(ax, 0.66, "WM:cortex", wm_ctx, 2.0, "#e45756", fmt_value(wm_ctx, 2), threshold=1.15)
    add_readout_bar(ax, 0.45, "Sub:cortex", sub_ctx, 2.0, "#f28e2b", fmt_value(sub_ctx, 2), threshold=1.25)
    add_readout_bar(ax, 0.24, "Motion max r", motion, 1.0, "#7f7f7f", fmt_value(motion, 2), threshold=0.35)
    ax.text(0.04, 0.07, peak_text, ha="left", va="bottom", fontsize=6.7, color="#333333")
    ax.text(0.96, 0.07, f"Cb:cortex {fmt_value(cereb_ctx, 2)}", ha="right", va="bottom", fontsize=6.7, color="#333333")


def movie_duration_seconds(movie: Path, ffprobe_bin: str) -> float:
    output = subprocess.check_output(
        [
            ffprobe_bin,
            "-v",
            "error",
            "-show_entries",
            "format=duration",
            "-of",
            "default=noprint_wrappers=1:nokey=1",
            str(movie),
        ],
        text=True,
    )
    duration = float(output.strip())
    if duration <= 0:
        raise RuntimeError(f"Movie duration must be positive for compression: {movie}")
    return duration


def compress_movie_to_target_size(
    source_movie: Path,
    output_movie: Path,
    target_size_mb: float,
    ffmpeg_bin: str,
    ffprobe_bin: str,
) -> None:
    duration = movie_duration_seconds(source_movie, ffprobe_bin)
    target_bits = float(target_size_mb) * 1024 * 1024 * 8
    video_kbps = max(50, int((target_bits / duration) * 0.97 / 1000))
    passlog = output_movie.with_suffix(".ffmpeg2pass")
    first_pass = [
        ffmpeg_bin,
        "-y",
        "-i",
        str(source_movie),
        "-an",
        "-c:v",
        "libx264",
        "-b:v",
        f"{video_kbps}k",
        "-pass",
        "1",
        "-passlogfile",
        str(passlog),
        "-f",
        "mp4",
        "/dev/null",
    ]
    second_pass = [
        ffmpeg_bin,
        "-y",
        "-i",
        str(source_movie),
        "-an",
        "-c:v",
        "libx264",
        "-b:v",
        f"{video_kbps}k",
        "-pass",
        "2",
        "-passlogfile",
        str(passlog),
        "-pix_fmt",
        "yuv420p",
        "-movflags",
        "+faststart",
        str(output_movie),
    ]
    subprocess.run(first_pass, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(second_pass, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for log_path in output_movie.parent.glob(passlog.name + "*"):
        try:
            log_path.unlink()
        except FileNotFoundError:
            pass


def frame_from_component(
    comp: np.ndarray,
    bg: np.ndarray | None,
    row: pd.Series,
    class_name: str,
    timecourse: np.ndarray | None,
    prior_labels: List[str],
    tr: float | None,
) -> np.ndarray:
    slice_idx = slice_indices_for_component(comp, n_slices=6)
    abs_comp = np.abs(comp)
    vmax = float(np.nanpercentile(abs_comp, 99.5)) if np.any(np.isfinite(abs_comp)) else 1.0
    vmax = max(vmax, np.finfo(float).eps)
    fig = plt.figure(figsize=(13.0, 7.2), dpi=120, facecolor="white")
    gs = fig.add_gridspec(4, 8, height_ratios=[1, 1, 1, 0.70], hspace=0.06, wspace=0.035)
    peak_text, psd_freq, psd_power, psd_peak = component_spectrum(timecourse, tr)
    title = (
        f"{component_title_class(row, prior_labels)} | component {int(row['component_id'])} | "
        f"k/r={float(row['log10_kappa_rho']):.2f} | NSI={float(row['decision_nsi']):.3f} | {peak_text}"
    )
    if np.isfinite(row.get("motion_max_abs_corr", np.nan)):
        title += f" | motion r={float(row['motion_max_abs_corr']):.2f}"
    fig.suptitle(title, fontsize=11, y=0.985)
    row_specs = [
        (0, "sag", "x", lambda i: (comp[i, :, :].T, None if bg is None else bg[i, :, :].T), slice_idx[0]),
        (1, "cor", "y", lambda i: (comp[:, i, :].T, None if bg is None else bg[:, i, :].T), slice_idx[1]),
        (2, "ax", "z", lambda i: (comp[:, :, i].T, None if bg is None else bg[:, :, i].T), slice_idx[2]),
    ]
    for grid_row, orient, dim_label, slicer, indices in row_specs:
        for grid_col, idx in enumerate(indices):
            ax = fig.add_subplot(gs[grid_row, grid_col])
            slc, bg_slc = slicer(int(idx))
            if bg_slc is not None:
                ax.imshow(bg_slc, cmap="gray", origin="lower", interpolation="nearest")
            else:
                ax.imshow(np.zeros_like(slc), cmap="gray", origin="lower", interpolation="nearest", vmin=0, vmax=1)
            masked = np.ma.masked_where(np.abs(slc) < 0.10 * vmax, slc)
            ax.imshow(masked, cmap="coolwarm", origin="lower", interpolation="nearest", vmin=-vmax, vmax=vmax, alpha=0.55)
            ax.text(
                0.02,
                0.95,
                f"{orient} {dim_label}={idx}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=6,
                color="white",
                bbox={"facecolor": "black", "alpha": 0.45, "pad": 1.5, "linewidth": 0},
            )
            ax.set_axis_off()
    add_decision_readout(fig.add_subplot(gs[0, 6:]), row)
    add_model_fit_readout(fig.add_subplot(gs[1, 6]), row, "rho")
    add_model_fit_readout(fig.add_subplot(gs[1, 7]), row, "kappa")
    add_spatial_context_readout(fig.add_subplot(gs[2, 6:]), row, peak_text)
    tc_ax = fig.add_subplot(gs[3, :6])
    if timecourse is not None and timecourse.size:
        tc = np.asarray(timecourse, dtype=float)
        tc = tc - np.nanmean(tc)
        sd = np.nanstd(tc)
        if sd > 0:
            tc = tc / sd
        tc_ax.plot(np.arange(tc.size), tc, color="#222222", linewidth=0.8)
        tc_ax.axhline(0, color="#999999", linewidth=0.6)
        tc_ax.set_ylabel("ICA TC (z)", fontsize=8)
        tc_ax.set_xlabel("Frame", fontsize=8)
        y = float(np.nanmax(np.abs(tc))) if np.any(np.isfinite(tc)) else 1.0
        tc_ax.set_ylim(-max(y, 1.0) * 1.08, max(y, 1.0) * 1.08)
    else:
        tc_ax.text(0.5, 0.5, "component time course unavailable", ha="center", va="center", fontsize=9)
        tc_ax.set_axis_off()
    tc_ax.tick_params(labelsize=7, length=2)
    tc_ax.grid(color="#dddddd", linewidth=0.4)
    psd_ax = fig.add_subplot(gs[3, 6:])
    if psd_freq is not None and psd_power is not None:
        mask = (psd_freq > 0) & (psd_freq <= 0.20)
        if np.any(mask):
            psd_ax.plot(psd_freq[mask], psd_power[mask], color="#333333", linewidth=0.8)
            if psd_peak is not None:
                psd_ax.axvline(psd_peak, color="#b2182b", linewidth=0.8)
            psd_ax.set_title(peak_text, fontsize=8)
            psd_ax.set_xlabel("Hz", fontsize=8)
            psd_ax.set_ylabel("Power", fontsize=8)
        else:
            psd_ax.text(0.5, 0.5, "spectrum unavailable", ha="center", va="center", fontsize=8)
    else:
        psd_ax.text(0.5, 0.5, "spectrum unavailable", ha="center", va="center", fontsize=8)
    psd_ax.tick_params(labelsize=7, length=2)
    psd_ax.grid(color="#dddddd", linewidth=0.4)
    fig.subplots_adjust(left=0.01, right=0.995, top=0.94, bottom=0.08)
    fig.canvas.draw()
    frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    plt.close(fig)
    return frame


def write_component_movie(
    class_name: str,
    ids: List[int],
    out: pd.DataFrame,
    comp_data: np.ndarray,
    bg: np.ndarray | None,
    timecourses: np.ndarray | None,
    prior_labels: List[str],
    tr: float | None,
    out_dir: Path,
    fps: float,
    write_mp4: bool,
    write_gif: bool,
    target_size_mb: float,
    keep_source_movie: bool,
    max_components: int,
    ffmpeg_bin: str,
    ffprobe_bin: str,
) -> Dict[str, object]:
    import io
    from PIL import Image

    out_dir.mkdir(parents=True, exist_ok=True)
    selected = [i for i in ids if 0 <= i < comp_data.shape[3]]
    if max_components > 0:
        selected = selected[:max_components]
    if not selected or (not write_mp4 and not write_gif):
        return {"class": class_name, "n_components": 0, "gif": "", "mp4": ""}
    gif_path = out_dir / f"{class_name}.gif"
    mp4_path = out_dir / f"{class_name}.mp4"
    source_mp4_path = out_dir / f"{class_name}.source.mp4"
    encode_mp4_path = source_mp4_path if target_size_mb > 0 else mp4_path
    wrote_mp4 = False
    wrote_gif = False
    mp4_proc = None
    gif_writer = None
    if not write_gif and gif_path.exists():
        gif_path.unlink()
    if write_mp4 and target_size_mb > 0 and not keep_source_movie and source_mp4_path.exists():
        source_mp4_path.unlink()
    if write_mp4 and shutil.which(ffmpeg_bin):
        mp4_proc = subprocess.Popen(
            [
                ffmpeg_bin,
                "-y",
                "-f",
                "image2pipe",
                "-vcodec",
                "png",
                "-r",
                str(float(fps)),
                "-i",
                "-",
                "-an",
                "-vcodec",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                str(encode_mp4_path),
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    if write_gif:
        import imageio.v2 as imageio

        gif_writer = imageio.get_writer(gif_path, mode="I", duration=1.0 / max(fps, np.finfo(float).eps))
    try:
        for cid in selected:
            frame = frame_from_component(
                comp_data[:, :, :, cid],
                bg,
                out.loc[out["component_id"] == cid].iloc[0],
                class_name,
                None if timecourses is None or cid >= timecourses.shape[1] else timecourses[:, cid],
                prior_labels,
                tr,
            )
            if gif_writer is not None:
                gif_writer.append_data(frame)
            if mp4_proc is not None and mp4_proc.stdin is not None:
                bio = io.BytesIO()
                Image.fromarray(frame).save(bio, format="PNG")
                mp4_proc.stdin.write(bio.getvalue())
    finally:
        if gif_writer is not None:
            gif_writer.close()
            wrote_gif = gif_path.is_file()
        if mp4_proc is not None:
            if mp4_proc.stdin is not None:
                mp4_proc.stdin.close()
            wrote_mp4 = mp4_proc.wait() == 0 and encode_mp4_path.is_file()
            if wrote_mp4 and target_size_mb > 0:
                try:
                    compress_movie_to_target_size(
                        source_mp4_path,
                        mp4_path,
                        target_size_mb,
                        ffmpeg_bin,
                        ffprobe_bin,
                    )
                    wrote_mp4 = mp4_path.is_file()
                except Exception:
                    wrote_mp4 = False
                    if mp4_path.exists():
                        mp4_path.unlink()
                    raise
                finally:
                    if not keep_source_movie and source_mp4_path.exists():
                        source_mp4_path.unlink()
    if write_mp4:
        if not wrote_mp4:
            if mp4_path.exists():
                mp4_path.unlink()
    return {
        "class": class_name,
        "n_components": len(selected),
        "gif": str(gif_path) if wrote_gif else "",
        "mp4": str(mp4_path) if wrote_mp4 else "",
    }


def write_qc_index(qc_dir: Path, records: List[Dict[str, object]]) -> None:
    rows = []
    for rec in records:
        gif = Path(str(rec.get("gif", ""))).name if rec.get("gif") else ""
        mp4 = Path(str(rec.get("mp4", ""))).name if rec.get("mp4") else ""
        media = ""
        if mp4:
            media = f'<video controls loop src="movies/{mp4}"></video>'
        elif gif:
            media = f'<img src="movies/{gif}" alt="{rec["class"]}">'
        rows.append(
            f"<section><h2>{rec['class']} ({rec['n_components']})</h2>{media}</section>"
        )
    html = """<!doctype html>
<html><head><meta charset="utf-8"><title>ME-ICA component reclassification QC</title>
<style>body{font-family:Arial,sans-serif;margin:24px;background:#f7f7f5;color:#222}section{margin:0 0 28px}video,img{max-width:100%;border:1px solid #ccc;background:#000}h1,h2{font-weight:600}</style>
</head><body><h1>ME-ICA component reclassification QC</h1>
""" + "\n".join(rows) + "\n</body></html>\n"
    (qc_dir / "index.html").write_text(html, encoding="utf-8")


def write_qc_class_outputs(out: pd.DataFrame, args: argparse.Namespace) -> None:
    if args.no_qc_class_outputs:
        return
    vol_path = args.component_nii if args.component_nii else args.tedana_dir / "betas_OC.nii.gz"
    if not vol_path.is_file():
        return
    qc_dir = args.out_dir / "ComponentClassQC"
    qc_dir.mkdir(parents=True, exist_ok=True)
    movie_dir = qc_dir / "movies"
    nifti_img = nib.load(str(vol_path))
    comp_data = np.asarray(nifti_img.dataobj, dtype=np.float32)
    cifti_path = args.tedana_dir / "betas_OC.dtseries.nii"
    bg = resample_bg_to_component(args.anat_bg, nifti_img)
    timecourses = load_component_timecourses(args.tedana_dir, comp_data.shape[3])
    prior_labels = load_prior_network_labels(args.priors_mat)
    tr = args.tr
    if tr is None and args.tr_file is not None and args.tr_file.is_file():
        vals = args.tr_file.read_text().split()
        if vals:
            tr = float(vals[0])
    records: List[Dict[str, object]] = []
    movie_classes = {
        item.strip()
        for item in str(args.qc_movie_classes).split(",")
        if item.strip()
    }
    for class_name, mask in class_masks(out).items():
        ids = out.loc[mask, "component_id"].astype(int).tolist()
        write_qc_ids(qc_dir / f"{class_name}_component_ids.txt", ids)
        write_nifti_subset(nifti_img, ids, qc_dir / f"betas_OC_{class_name}.nii.gz")
        write_cifti_subset(cifti_path, ids, qc_dir / f"betas_OC_{class_name}.dscalar.nii")
        if args.qc_movies and (not movie_classes or class_name in movie_classes):
            records.append(
                write_component_movie(
                    class_name,
                    ids,
                    out,
                    comp_data,
                    bg,
                    timecourses,
                    prior_labels,
                    tr,
                    movie_dir,
                    fps=args.qc_movie_fps,
                    write_mp4=bool(args.qc_movie_mp4),
                    write_gif=bool(args.qc_movie_gif),
                    target_size_mb=float(args.qc_movie_target_size_mb),
                    keep_source_movie=bool(args.qc_movie_keep_source_movie),
                    max_components=args.qc_movie_max_components,
                    ffmpeg_bin=args.ffmpeg_bin,
                    ffprobe_bin=args.ffprobe_bin,
                )
            )
    if records:
        write_qc_index(qc_dir, records)
        (qc_dir / "qc_media_summary.json").write_text(json.dumps(records, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tedana-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--training-features", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--priors-mat", type=Path, default=DEFAULT_PRIORS)
    parser.add_argument("--neighbor-table", type=Path, default=DEFAULT_NEIGHBORS)
    parser.add_argument("--component-nii", type=Path, default=None)
    parser.add_argument("--anat-bg", type=Path, default=None)
    parser.add_argument("--tr", type=float, default=None)
    parser.add_argument("--tr-file", type=Path, default=None)
    parser.add_argument("--tc-peak-max-hz", type=float, default=0.20)
    parser.add_argument("--motion-params", type=Path, default=None)
    parser.add_argument("--motion-corr-threshold", type=float, default=0.35)
    parser.add_argument("--aseg", type=Path, default=None)
    parser.add_argument("--wm-mask", type=Path, default=None)
    parser.add_argument("--cortical-ribbon", type=Path, default=None)
    parser.add_argument("--subcortical-rois", type=Path, default=None)
    parser.add_argument("--ridge-lambda", type=float, default=10.0)
    parser.add_argument("--spectral-kmax", type=int, default=300)
    parser.add_argument("--rf-trees", type=int, default=200)
    parser.add_argument("--rf-min-samples-leaf", type=int, default=5)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--kill-mode", choices=("adaptive", "fixed"), default="adaptive")
    parser.add_argument("--kill-nsi", type=float, default=0.05)
    parser.add_argument("--kill-nsi-min", type=float, default=0.04)
    parser.add_argument("--kill-nsi-max", type=float, default=0.10)
    parser.add_argument("--kill-intercept", type=float, default=0.14)
    parser.add_argument("--kill-slope", type=float, default=0.25)
    parser.add_argument("--decision-nsi-col", choices=("nsi", "nsi_fc_demeaned"), default="nsi_fc_demeaned")
    parser.add_argument("--pred-kill-logratio", type=float, default=0.15)
    parser.add_argument("--pred-protect-logratio", type=float, default=0.45)
    parser.add_argument("--observed-protect-logratio", type=float, default=0.50)
    parser.add_argument("--protect-template-corr", type=float, default=0.18)
    parser.add_argument("--protect-morans", type=float, default=0.85)
    parser.add_argument("--protect-spectral-slope-max", type=float, default=-0.45)
    parser.add_argument("--protect-demeaned-nsi", type=float, default=0.10)
    parser.add_argument("--wm-kill-ratio", type=float, default=1.15)
    parser.add_argument("--subcort-protect-ratio", type=float, default=1.25)
    parser.add_argument("--cerebellum-protect-ratio", type=float, default=1.25)
    parser.add_argument("--rescue-nsi", type=float, default=0.20)
    parser.add_argument("--pred-rescue-logratio", type=float, default=0.35)
    parser.add_argument("--guardrail-kappa-rho", type=int, default=1)
    parser.add_argument("--kill-priority-enable", type=int, default=1)
    parser.add_argument("--kill-priority-w-nsi", type=float, default=0.35)
    parser.add_argument("--kill-priority-w-pred", type=float, default=0.45)
    parser.add_argument("--kill-priority-w-var", type=float, default=0.20)
    parser.add_argument("--kill-max-count", type=int, default=0)
    parser.add_argument("--no-qc-class-outputs", action="store_true")
    parser.add_argument("--qc-movies", type=int, default=0)
    parser.add_argument("--qc-movie-classes", default="")
    parser.add_argument("--qc-movie-max-components", type=int, default=80)
    parser.add_argument("--qc-movie-fps", type=float, default=0.35)
    parser.add_argument("--qc-movie-mp4", type=int, default=1)
    parser.add_argument("--qc-movie-gif", type=int, default=0)
    parser.add_argument("--qc-movie-target-size-mb", type=float, default=0.0)
    parser.add_argument("--qc-movie-keep-source-movie", type=int, default=0)
    parser.add_argument("--ffmpeg-bin", default="ffmpeg")
    parser.add_argument("--ffprobe-bin", default="ffprobe")
    parser.add_argument("--copy-key-files", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    model = train_spatial_model(args.training_features, args)
    df = compute_run_features(args)
    out = refine(df, model, args)
    out = add_timecourse_spectral_features(out, args)
    out = add_motion_coupling_features(out, args)
    out.to_csv(args.out_dir / "ComponentDecisions.tsv", sep="\t", index=False)
    accepted = out.loc[out["accepted_final"], "component_id"].astype(int).tolist()
    rejected = out.loc[~out["accepted_final"], "component_id"].astype(int).tolist()
    write_ids(args.out_dir / "AcceptedComponents.txt", accepted)
    write_ids(args.out_dir / "RejectedComponents.txt", rejected)
    summary = summarize(out, args)
    (args.out_dir / "ClassificationSummary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.out_dir / "ClassificationSummary.txt").write_text(
        "\n".join(f"{k}: {v}" for k, v in summary.items() if k != "settings") + "\n",
        encoding="utf-8",
    )
    make_traffic_plot(out, args.out_dir / "KappaRhoLog_vs_NSI_spatial_traffic.png", args)
    make_predicted_observed_plot(out, args.out_dir / "Predicted_vs_Observed_LogKappaRho.png")
    write_qc_class_outputs(out, args)
    if args.copy_key_files:
        for name in ["betas_OC.dtseries.nii", "desc-ICA_decomposition.json", "desc-tedana_metrics.tsv", "ica_decomposition.json", "ica_metrics.tsv"]:
            src = args.tedana_dir / name
            if src.exists():
                shutil.copy2(src, args.out_dir / name)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

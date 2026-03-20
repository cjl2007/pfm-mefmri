#!/usr/bin/env python3
"""ME-ICA reclassification and manual-override parser.

Outputs tedana-compatible component lists:
- AcceptedComponents.txt
- RejectedComponents.txt
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from nibabel.cifti2 import cifti2_axes


def _component_id(component_value: str) -> int:
    m = re.search(r"(\d+)$", str(component_value))
    if not m:
        raise ValueError(f"Could not parse component id from: {component_value}")
    return int(m.group(1))


def _load_component_metrics_table(ted_dir: Path) -> pd.DataFrame:
    """Load tedana component metrics across v12/v10 output layouts."""
    for metrics_tsv in (ted_dir / "ica_metrics.tsv", ted_dir / "desc-tedana_metrics.tsv"):
        if not metrics_tsv.is_file():
            continue
        df = pd.read_csv(metrics_tsv, sep="\t")
        for c in ("Component", "kappa", "rho", "classification"):
            if c not in df.columns:
                raise ValueError(f"Missing required column '{c}' in {metrics_tsv}")
        df["component_id"] = df["Component"].map(_component_id).astype(int)
        return df.sort_values("component_id").reset_index(drop=True)

    # Some versions emit comp_table.tsv naming.
    comp_table_tsv = ted_dir / "comp_table.tsv"
    if comp_table_tsv.is_file():
        df = pd.read_csv(comp_table_tsv, sep="\t")
        rename_map = {}
        if "component" in df.columns and "Component" not in df.columns:
            rename_map["component"] = "Component"
        if "Kappa" in df.columns and "kappa" not in df.columns:
            rename_map["Kappa"] = "kappa"
        if "Rho" in df.columns and "rho" not in df.columns:
            rename_map["Rho"] = "rho"
        if "Classification" in df.columns and "classification" not in df.columns:
            rename_map["Classification"] = "classification"
        if rename_map:
            df = df.rename(columns=rename_map)
        for c in ("Component", "kappa", "rho", "classification"):
            if c not in df.columns:
                raise ValueError(f"Missing required column '{c}' in {comp_table_tsv}")
        df["component_id"] = df["Component"].map(_component_id).astype(int)
        return df.sort_values("component_id").reset_index(drop=True)

    # v10 fallback: parse per-component entries from ica_decomposition.json
    for decomp_json in (ted_dir / "ica_decomposition.json", ted_dir / "desc-ICA_decomposition.json"):
        if not decomp_json.is_file():
            continue
        obj = json.loads(decomp_json.read_text())
        rows: list[dict] = []
        for key, val in obj.items():
            if not str(key).startswith("ica_"):
                continue
            try:
                cid = _component_id(str(key))
            except ValueError:
                continue
            if not isinstance(val, dict):
                continue
            row = {"Component": str(key), "component_id": int(cid)}
            for metric in ("classification", "kappa", "rho", "variance explained", "normalized variance explained"):
                if metric in val:
                    row[metric] = val[metric]
            rows.append(row)
        if not rows:
            raise ValueError(f"No ICA component rows found in {decomp_json}")
        df = pd.DataFrame(rows)
        for c in ("classification", "kappa", "rho"):
            if c not in df.columns:
                raise ValueError(f"Missing required key '{c}' in {decomp_json}")
        return df.sort_values("component_id").reset_index(drop=True)

    raise FileNotFoundError(
        f"Missing tedana component table in {ted_dir} "
        f"(tried ica_metrics.tsv, desc-tedana_metrics.tsv, comp_table.tsv, "
        f"ica_decomposition.json, desc-ICA_decomposition.json)"
    )


def _read_manual_component_ids(fig_dir: Path) -> set[int]:
    if not fig_dir.is_dir():
        return set()
    ids: set[int] = set()
    for p in fig_dir.glob("*.png"):
        m = re.search(r"(?:comp|ica)[_\-]?(\d+)", p.name, flags=re.IGNORECASE)
        if not m:
            continue
        ids.add(int(m.group(1)))
    return ids


def _mat_struct_get(s: np.ndarray, field: str):
    if isinstance(s, np.ndarray) and s.dtype.names and field in s.dtype.names:
        return s[field][0, 0]
    raise ValueError(f"Missing field '{field}' in Priors struct")


def _load_priors_fc(priors_mat: Path) -> np.ndarray:
    m = sio.loadmat(priors_mat)
    if "Priors" not in m:
        raise ValueError(f"{priors_mat} missing variable 'Priors'")
    pri = m["Priors"]
    fc = _mat_struct_get(pri, "FC")
    return np.asarray(fc, dtype=np.float64)


def _ridge_r2(maps: np.ndarray, priors_fc: np.ndarray, lam: float = 10.0) -> np.ndarray:
    # maps: n_vertices x n_components
    # priors_fc: n_vertices x n_templates
    u, s, vt = np.linalg.svd(priors_fc, full_matrices=False)
    v = vt.T
    ut_y = u.T @ maps
    w = s / (s**2 + lam)
    b = v @ (w[:, None] * ut_y)
    yhat = priors_fc @ b
    sse = np.sum((maps - yhat) ** 2, axis=0)
    sst = np.sum((maps - np.mean(maps, axis=0, keepdims=True)) ** 2, axis=0)
    r2 = 1 - sse / np.maximum(sst, np.finfo(np.float64).eps)
    r2[r2 < 0] = np.nan
    return r2


def _compute_nsi_from_cifti(betas_cifti: Path, priors_mat: Path, n_components: int) -> np.ndarray:
    priors_fc = _load_priors_fc(priors_mat)
    img = nib.load(str(betas_cifti))
    data = np.asarray(img.dataobj, dtype=np.float64)

    # Expect either (grayordinates, components) or (components, grayordinates).
    if data.ndim != 2:
        raise ValueError(f"Expected 2D CIFTI data, got shape {data.shape}")
    if data.shape[1] == n_components:
        maps = data
    elif data.shape[0] == n_components:
        maps = data.T
    else:
        raise ValueError(
            f"CIFTI shape {data.shape} does not match component count {n_components}"
        )

    if maps.shape[0] < priors_fc.shape[0]:
        raise ValueError(
            f"CIFTI rows ({maps.shape[0]}) < Priors.FC rows ({priors_fc.shape[0]})"
        )

    # Legacy convention uses cortical rows aligned to Priors.FC.
    cortical_maps = maps[: priors_fc.shape[0], :]
    nsi = _ridge_r2(cortical_maps, priors_fc, lam=10.0)
    nsi = np.where(np.isnan(nsi), 0.0, nsi)
    return nsi


def _compute_template_rho_from_cifti(betas_cifti: Path, priors_mat: Path, n_components: int) -> np.ndarray:
    """Legacy-style max abs corr with network templates (Priors.FC)."""
    priors_fc = _load_priors_fc(priors_mat)
    img = nib.load(str(betas_cifti))
    data = np.asarray(img.dataobj, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D CIFTI data, got shape {data.shape}")
    if data.shape[1] == n_components:
        maps = data
    elif data.shape[0] == n_components:
        maps = data.T
    else:
        raise ValueError(
            f"CIFTI shape {data.shape} does not match component count {n_components}"
        )
    if maps.shape[0] < priors_fc.shape[0]:
        raise ValueError(
            f"CIFTI rows ({maps.shape[0]}) < Priors.FC rows ({priors_fc.shape[0]})"
        )
    x = maps[: priors_fc.shape[0], :]
    # Correlation between each component map and each template.
    x0 = x - np.nanmean(x, axis=0, keepdims=True)
    t0 = priors_fc - np.nanmean(priors_fc, axis=0, keepdims=True)
    num = t0.T @ x0
    den = np.sqrt(
        np.sum(t0 * t0, axis=0, keepdims=True).T
        * np.sum(x0 * x0, axis=0, keepdims=True)
    )
    r = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return np.max(np.abs(r), axis=0)


def _compute_subcortical_ratio_from_cifti(betas_cifti: Path, n_components: int) -> np.ndarray:
    """Compute max subcortical:cortical mean(|loading|) ratio per component."""
    img = nib.load(str(betas_cifti))
    data = np.asarray(img.dataobj, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D CIFTI data, got shape {data.shape}")

    if data.shape[1] == n_components:
        maps = data
    elif data.shape[0] == n_components:
        maps = data.T
    else:
        raise ValueError(
            f"CIFTI shape {data.shape} does not match component count {n_components}"
        )

    ax0 = img.header.get_axis(0)
    ax1 = img.header.get_axis(1)
    bm_axis = None
    if isinstance(ax0, cifti2_axes.BrainModelAxis) and len(ax0) == maps.shape[0]:
        bm_axis = ax0
    elif isinstance(ax1, cifti2_axes.BrainModelAxis) and len(ax1) == maps.shape[0]:
        bm_axis = ax1
    if bm_axis is None:
        raise ValueError("Could not align BrainModelAxis with CIFTI grayordinate dimension")

    cortex_idx = np.zeros(maps.shape[0], dtype=bool)
    subcort_masks = {
        "accumbens": np.zeros(maps.shape[0], dtype=bool),
        "caudate": np.zeros(maps.shape[0], dtype=bool),
        "pallidum": np.zeros(maps.shape[0], dtype=bool),
        "putamen": np.zeros(maps.shape[0], dtype=bool),
        "thalamus": np.zeros(maps.shape[0], dtype=bool),
        "hippocampus": np.zeros(maps.shape[0], dtype=bool),
        "amygdala": np.zeros(maps.shape[0], dtype=bool),
        "cerebellum": np.zeros(maps.shape[0], dtype=bool),
    }
    subcort_names = {f"cifti_structure_{k}" for k in subcort_masks.keys()}
    subcort_prefixes = tuple(subcort_names)

    for struct_name, slc, _ in bm_axis.iter_structures():
        s = struct_name.lower()
        if s in ("cifti_structure_cortex_left", "cifti_structure_cortex_right"):
            cortex_idx[slc] = True
            continue
        if not s.startswith(subcort_prefixes):
            continue
        for key in subcort_masks:
            if s.startswith(f"cifti_structure_{key}"):
                subcort_masks[key][slc] = True
                break

    n_components_eff = maps.shape[1]
    out = np.full(n_components_eff, np.nan, dtype=np.float64)
    if not np.any(cortex_idx):
        return out

    for i in range(n_components_eff):
        comp = np.abs(maps[:, i])
        cort_mean = float(np.mean(comp[cortex_idx]))
        if cort_mean <= 0:
            out[i] = np.inf
            continue
        ratios = []
        for _, mask in subcort_masks.items():
            if np.any(mask):
                ratios.append(float(np.mean(comp[mask]) / (cort_mean + np.finfo(np.float64).eps)))
        if ratios:
            out[i] = float(np.max(ratios))
    return out


def _initial_keep_from_classification(df: pd.DataFrame) -> pd.Series:
    """Map tedana class labels to initial keep/reject across versions.

    Legacy tedana v0.x often includes `ignored`; these should not be treated
    as initially accepted.
    """
    cls = df["classification"].astype(str).str.strip().str.lower()
    keep_labels = {"accepted", "provisionalaccept"}
    reject_labels = {"rejected", "ignored", "provisionalreject"}
    keep = cls.isin(keep_labels)
    unknown = ~(cls.isin(keep_labels) | cls.isin(reject_labels))
    if unknown.any():
        # Conservative default for unknown labels.
        keep.loc[unknown] = False
    return keep


def _nsi_refine(
    df: pd.DataFrame,
    nsi: np.ndarray,
    subcort_ratio: np.ndarray | None,
    *,
    guardrail: bool,
    rescue_nsi: float,
    rescue_quantile: float,
    kill_mode: str,
    kill_nsi: float,
    kill_nsi_min: float,
    kill_nsi_max: float,
    kill_intercept: float,
    kill_slope: float,
    subcort_ratio_thresh: float,
    kill_priority_enable: bool,
    kill_priority_w_logratio: float,
    kill_priority_w_nsi: float,
    kill_priority_w_var: float,
    kill_var_floor_quantile: float,
    kill_cumvar_cap: float,
) -> pd.DataFrame:
    out = df.copy()
    out["NSI"] = nsi
    keep_init = _initial_keep_from_classification(out)
    if guardrail:
        keep_init &= (out["kappa"] > out["rho"])
    out["keep_init"] = keep_init
    if subcort_ratio is not None and len(subcort_ratio) == len(out):
        out["subcort_ratio"] = subcort_ratio
    else:
        out["subcort_ratio"] = np.nan

    keep_final = keep_init.to_numpy(copy=True)
    rescued = np.zeros(len(out), dtype=bool)
    killed = np.zeros(len(out), dtype=bool)
    kill_candidate = np.zeros(len(out), dtype=bool)

    sig_nsi = out.loc[keep_init & np.isfinite(out["NSI"]), "NSI"].to_numpy(dtype=float)
    nsi_floor = np.quantile(sig_nsi, rescue_quantile) if sig_nsi.size else np.inf

    eps = np.finfo(np.float64).eps
    for i in range(len(out)):
        nsi_i = float(out.iloc[i]["NSI"])
        if not np.isfinite(nsi_i):
            continue
        kappa_i = float(out.iloc[i]["kappa"])
        rho_i = float(out.iloc[i]["rho"])

        want_rescue = (
            (not bool(keep_init.iloc[i]))
            and (nsi_i >= rescue_nsi)
            and (nsi_i >= nsi_floor)
        )
        if guardrail:
            want_rescue = want_rescue and (kappa_i > rho_i)

        if kill_mode == "adaptive":
            logratio = np.log10((kappa_i / max(rho_i, eps)) + eps)
            nsi_kill_thresh = np.clip(
                kill_intercept - kill_slope * logratio, kill_nsi_min, kill_nsi_max
            )
        else:
            nsi_kill_thresh = kill_nsi
        want_kill = bool(keep_init.iloc[i]) and (nsi_i < nsi_kill_thresh)

        # Subcortical guardrail (kill-only):
        # if tedana initially accepted and kappa > rho, do not kill
        # components with strongly subcortical loading.
        if want_kill and np.isfinite(out.iloc[i]["subcort_ratio"]) and (kappa_i > rho_i):
            if float(out.iloc[i]["subcort_ratio"]) >= subcort_ratio_thresh:
                want_kill = False

        if want_rescue:
            keep_final[i] = True
            rescued[i] = True
        elif want_kill:
            kill_candidate[i] = True

    # Optional kill prioritization: choose a principled subset of kill candidates.
    if kill_priority_enable and np.any(kill_candidate):
        idx = np.where(kill_candidate)[0]
        # Defensive fallback if variance is unavailable.
        if "variance explained" in out.columns:
            v = out["variance explained"].to_numpy(dtype=float)
        elif "normalized variance explained" in out.columns:
            v = out["normalized variance explained"].to_numpy(dtype=float)
        else:
            v = np.ones(len(out), dtype=float)
        v = np.where(np.isfinite(v), v, 0.0)

        # Evidence floor on candidate variance within kill set.
        cand_v = v[idx]
        floor = float(np.quantile(cand_v, kill_var_floor_quantile)) if len(cand_v) else 0.0
        pool_mask = cand_v >= floor
        pool_idx = idx[pool_mask]
        if len(pool_idx) == 0:
            pool_idx = idx

        # Score = lower logratio + lower NSI + higher variance.
        def z(a: np.ndarray) -> np.ndarray:
            a = np.asarray(a, dtype=float)
            s = float(np.std(a))
            if s <= 0:
                return np.zeros_like(a)
            return (a - float(np.mean(a))) / s

        logratio = np.log10(
            (out.loc[pool_idx, "kappa"].to_numpy(dtype=float) /
             np.maximum(out.loc[pool_idx, "rho"].to_numpy(dtype=float), eps)) + eps
        )
        nsi_pool = out.loc[pool_idx, "NSI"].to_numpy(dtype=float)
        var_pool = v[pool_idx]
        score = (
            kill_priority_w_logratio * z(-logratio) +
            kill_priority_w_nsi * z(-nsi_pool) +
            kill_priority_w_var * z(var_pool)
        )
        order = np.argsort(-score)
        ranked_idx = pool_idx[order]
        ranked_var = var_pool[order]

        # Cumulative variance cap inside prioritized pool.
        total_var = float(np.sum(ranked_var))
        if total_var <= 0:
            selected = ranked_idx
        else:
            target = kill_cumvar_cap * total_var
            cs = np.cumsum(ranked_var)
            n_keep = int(np.sum(cs <= target))
            if n_keep < len(ranked_idx):
                n_keep += 1
            selected = ranked_idx[:n_keep]

        killed[:] = False
        killed[selected] = True
    else:
        killed = kill_candidate.copy()

    keep_final[killed] = False

    if guardrail:
        keep_final &= (out["kappa"].to_numpy() > out["rho"].to_numpy())

    out["kill_candidate"] = kill_candidate
    out["rescued"] = rescued
    out["killed"] = killed
    out["accepted_final"] = keep_final
    return out


def _legacy_template_rho_refine(
    df: pd.DataFrame,
    template_rho: np.ndarray,
    rho_rescue: float,
    rho_reject: float,
) -> pd.DataFrame:
    out = df.copy()
    keep_init = _initial_keep_from_classification(out)
    rho_t = np.asarray(template_rho, dtype=float)
    if rho_t.shape[0] != len(out):
        raise ValueError("template_rho length mismatch")
    rescued = (~keep_init.to_numpy()) & (rho_t > rho_rescue)
    killed = (keep_init.to_numpy()) & (rho_t < rho_reject)
    keep_final = keep_init.to_numpy()
    keep_final[rescued] = True
    keep_final[killed] = False
    out["template_rho"] = rho_t
    out["keep_init"] = keep_init
    out["rescued"] = rescued
    out["killed"] = killed
    out["accepted_final"] = keep_final
    return out


def _write_ids(path: Path, ids: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(str(i) for i in ids) + ("\n" if ids else ""))

def _make_nsi_plot(
    out: pd.DataFrame,
    out_png: Path,
    *,
    kill_mode: str,
    kill_nsi: float,
    kill_nsi_min: float,
    kill_nsi_max: float,
    kill_intercept: float,
    kill_slope: float,
) -> None:
    if "NSI" not in out.columns:
        return
    eps = np.finfo(np.float64).eps
    x = np.log10((out["kappa"].to_numpy(dtype=float) / np.maximum(out["rho"].to_numpy(dtype=float), eps)) + eps)
    y = out["NSI"].to_numpy(dtype=float)

    keep_init = out["keep_init"].to_numpy(dtype=bool)
    rescued = out["rescued"].to_numpy(dtype=bool)
    killed = out["killed"].to_numpy(dtype=bool)

    plt.figure(figsize=(7, 4), dpi=150)
    plt.scatter(x[~keep_init], y[~keep_init], s=12, c="tab:red", alpha=0.7, label="init reject")
    plt.scatter(x[keep_init], y[keep_init], s=12, c="tab:green", alpha=0.7, label="init accept")
    if np.any(rescued):
        plt.scatter(x[rescued], y[rescued], s=42, c="gold", edgecolors="k", linewidths=0.4, label="rescued")
    if np.any(killed):
        plt.scatter(x[killed], y[killed], s=42, c="cyan", edgecolors="k", linewidths=0.4, label="killed")
    if np.any(np.isfinite(x)):
        x_line = np.linspace(np.nanmin(x), np.nanmax(x), 256)
        if kill_mode == "adaptive":
            y_line = np.clip(kill_intercept - kill_slope * x_line, kill_nsi_min, kill_nsi_max)
            plt.plot(x_line, y_line, "k--", linewidth=1.1, label="adaptive kill threshold")
        else:
            plt.axhline(kill_nsi, color="k", linestyle="--", linewidth=1.1, label="fixed kill threshold")
    plt.xlabel("log10(kappa / rho)")
    plt.ylabel("NSI")
    plt.title("NSI vs log(kappa/rho)")
    plt.grid(alpha=0.25)
    plt.legend(loc="best", frameon=False)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Reclassify tedana components.")
    ap.add_argument("--data-dir", required=True, help="Run directory containing Tedana/")
    ap.add_argument("--tedana-dir", default="", help="Explicit Tedana output directory (default: <data-dir>/Tedana)")
    ap.add_argument("--out-dir", default="", help="Output dir (default: <data-dir>/Tedana+ManualComponentClassification)")
    ap.add_argument(
        "--classifier-mode",
        default="nsi",
        choices=["nsi", "legacy_template_rho", "none"],
    )
    ap.add_argument("--priors-mat", default="", help="Priors.mat for NSI mode")
    ap.add_argument("--betas-cifti", default="", help="betas_OC.dtseries.nii for NSI mode")
    ap.add_argument("--rho-rescue", type=float, default=0.30)
    ap.add_argument("--rho-reject", type=float, default=0.10)
    ap.add_argument("--rescue-nsi", type=float, default=0.20)
    ap.add_argument("--rescue-quantile", type=float, default=0.10)
    ap.add_argument("--kill-mode", default="adaptive", choices=["adaptive", "fixed"])
    ap.add_argument("--kill-nsi", type=float, default=0.05)
    ap.add_argument("--kill-nsi-min", type=float, default=0.02)
    ap.add_argument("--kill-nsi-max", type=float, default=0.10)
    ap.add_argument("--kill-intercept", type=float, default=0.10)
    ap.add_argument("--kill-slope", type=float, default=0.50)
    ap.add_argument("--guardrail-kappa-rho", type=int, default=1)
    ap.add_argument("--subcort-ratio-thresh", type=float, default=5.0)
    ap.add_argument("--kill-priority-enable", type=int, default=1)
    ap.add_argument("--kill-priority-w-logratio", type=float, default=0.50)
    ap.add_argument("--kill-priority-w-nsi", type=float, default=0.30)
    ap.add_argument("--kill-priority-w-var", type=float, default=0.20)
    ap.add_argument("--kill-var-floor-quantile", type=float, default=0.60)
    ap.add_argument("--kill-cumvar-cap", type=float, default=0.95)
    args = ap.parse_args()

    run_dir = Path(args.data_dir).resolve()
    ted_dir = Path(args.tedana_dir).resolve() if args.tedana_dir else (run_dir / "Tedana")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (run_dir / "Tedana+ManualComponentClassification")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = _load_component_metrics_table(ted_dir)

    man_rej = _read_manual_component_ids(ted_dir / "figures" / "ManuallyRejected")
    man_acc = _read_manual_component_ids(ted_dir / "figures" / "ManuallyAccepted")

    mode_used = args.classifier_mode
    if args.classifier_mode == "nsi":
        betas_cifti = Path(args.betas_cifti) if args.betas_cifti else (ted_dir / "betas_OC.dtseries.nii")
        priors_mat = Path(args.priors_mat) if args.priors_mat else None
        if priors_mat and priors_mat.is_file() and betas_cifti.is_file():
            nsi = _compute_nsi_from_cifti(betas_cifti, priors_mat, len(df))
            subcort_ratio = _compute_subcortical_ratio_from_cifti(betas_cifti, len(df))
            out = _nsi_refine(
                df,
                nsi,
                subcort_ratio,
                guardrail=bool(args.guardrail_kappa_rho),
                rescue_nsi=args.rescue_nsi,
                rescue_quantile=args.rescue_quantile,
                kill_mode=args.kill_mode,
                kill_nsi=args.kill_nsi,
                kill_nsi_min=args.kill_nsi_min,
                kill_nsi_max=args.kill_nsi_max,
                kill_intercept=args.kill_intercept,
                kill_slope=args.kill_slope,
                subcort_ratio_thresh=args.subcort_ratio_thresh,
                kill_priority_enable=bool(args.kill_priority_enable),
                kill_priority_w_logratio=args.kill_priority_w_logratio,
                kill_priority_w_nsi=args.kill_priority_w_nsi,
                kill_priority_w_var=args.kill_priority_w_var,
                kill_var_floor_quantile=args.kill_var_floor_quantile,
                kill_cumvar_cap=args.kill_cumvar_cap,
            )
        else:
            mode_used = "none_fallback"
            out = df.copy()
            keep_init = _initial_keep_from_classification(out)
            out["keep_init"] = keep_init
            out["rescued"] = False
            out["killed"] = False
            out["accepted_final"] = keep_init.to_numpy()
    elif args.classifier_mode == "legacy_template_rho":
        betas_cifti = Path(args.betas_cifti) if args.betas_cifti else (ted_dir / "betas_OC.dtseries.nii")
        priors_mat = Path(args.priors_mat) if args.priors_mat else None
        if priors_mat and priors_mat.is_file() and betas_cifti.is_file():
            trho = _compute_template_rho_from_cifti(betas_cifti, priors_mat, len(df))
            out = _legacy_template_rho_refine(df, trho, args.rho_rescue, args.rho_reject)
        else:
            mode_used = "none_fallback"
            out = df.copy()
            keep_init = _initial_keep_from_classification(out)
            out["keep_init"] = keep_init
            out["rescued"] = False
            out["killed"] = False
            out["accepted_final"] = keep_init.to_numpy()
    else:
        out = df.copy()
        keep_init = _initial_keep_from_classification(out)
        out["keep_init"] = keep_init
        out["rescued"] = False
        out["killed"] = False
        out["accepted_final"] = keep_init.to_numpy()

    out["manual_accept"] = out["component_id"].isin(man_acc)
    out["manual_reject"] = out["component_id"].isin(man_rej)

    # Manual overrides take precedence.
    accepted_final = out["accepted_final"].to_numpy(copy=True)
    accepted_final[out["manual_accept"].to_numpy()] = True
    accepted_final[out["manual_reject"].to_numpy()] = False
    out["accepted_final"] = accepted_final

    rescued_count = int(np.sum(out["rescued"].to_numpy(dtype=bool)))
    killed_count = int(np.sum(out["killed"].to_numpy(dtype=bool)))
    kill_candidate_count = int(np.sum(out.get("kill_candidate", pd.Series(False, index=out.index)).to_numpy(dtype=bool)))
    init_accept_count = int(np.sum(out["keep_init"].to_numpy(dtype=bool)))
    init_reject_count = int(len(out) - init_accept_count)
    changed_count = int(np.sum(out["keep_init"].to_numpy(dtype=bool) != out["accepted_final"].to_numpy(dtype=bool)))

    accepted_ids = out.loc[out["accepted_final"], "component_id"].astype(int).tolist()
    rejected_ids = out.loc[~out["accepted_final"], "component_id"].astype(int).tolist()
    accepted_ids = sorted(set(accepted_ids))
    rejected_ids = sorted(set(rejected_ids))

    _write_ids(out_dir / "AcceptedComponents.txt", accepted_ids)
    _write_ids(out_dir / "RejectedComponents.txt", rejected_ids)
    out.to_csv(out_dir / "ComponentDecisions.tsv", sep="\t", index=False)

    summary = {
        "n_components": int(len(out)),
        "mode_requested": args.classifier_mode,
        "mode_used": mode_used,
        "manual_accept_count": len(man_acc),
        "manual_reject_count": len(man_rej),
        "initial_accept_count": init_accept_count,
        "initial_reject_count": init_reject_count,
        "rescued_count": rescued_count,
        "killed_count": killed_count,
        "kill_candidate_count": kill_candidate_count,
        "changed_count": changed_count,
        "accepted_count": len(accepted_ids),
        "rejected_count": len(rejected_ids),
    }
    (out_dir / "ClassificationSummary.json").write_text(json.dumps(summary, indent=2))
    summary_txt = [
        f"mode_requested: {summary['mode_requested']}",
        f"mode_used: {summary['mode_used']}",
        f"n_components: {summary['n_components']}",
        f"initial_accept_count: {summary['initial_accept_count']}",
        f"initial_reject_count: {summary['initial_reject_count']}",
        f"rescued_count: {summary['rescued_count']}",
        f"killed_count: {summary['killed_count']}",
        f"kill_candidate_count: {summary['kill_candidate_count']}",
        f"changed_count: {summary['changed_count']}",
        f"manual_accept_count: {summary['manual_accept_count']}",
        f"manual_reject_count: {summary['manual_reject_count']}",
        f"accepted_count: {summary['accepted_count']}",
        f"rejected_count: {summary['rejected_count']}",
    ]
    (out_dir / "ClassificationSummary.txt").write_text("\n".join(summary_txt) + "\n")

    if mode_used.startswith("nsi") and "NSI" in out.columns:
        _make_nsi_plot(
            out,
            out_dir / "KappaRhoLog_vs_NSI.png",
            kill_mode=args.kill_mode,
            kill_nsi=args.kill_nsi,
            kill_nsi_min=args.kill_nsi_min,
            kill_nsi_max=args.kill_nsi_max,
            kill_intercept=args.kill_intercept,
            kill_slope=args.kill_slope,
        )

    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Pilot model: predict tedana kappa/rho metrics from ICA-map NSI.

This script is intentionally self-contained and writes only to the current
project by default. It reads TMS-PFM Tedana outputs and ME-fMRI pipeline priors.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import nibabel as nib
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler


DEFAULT_PRIORS = Path(__file__).resolve().parents[1] / "res0urces" / "Priors.mat"
DEFAULT_NEIGHBORS = (
    Path(__file__).resolve().parents[1] / "lib" / "pfm-nsi" /
    "pfm_nsi" / "models" / "cifti_surf_neighbors_lr_normalwall.npz"
)


def component_id(component_value: str) -> int:
    match = re.search(r"(\d+)$", str(component_value))
    if not match:
        raise ValueError("Could not parse component id from {}".format(component_value))
    return int(match.group(1))


def scalar_or_nan(value) -> float:
    try:
        if isinstance(value, (str, bytes)):
            return float(value)
        if np.isscalar(value):
            return float(value)
    except Exception:
        return float("nan")
    return float("nan")


def mat_struct_get(s: np.ndarray, field: str):
    if isinstance(s, np.ndarray) and s.dtype.names and field in s.dtype.names:
        return s[field][0, 0]
    raise ValueError("Missing field '{}' in Priors struct".format(field))


def load_priors_fc(priors_mat: Path) -> np.ndarray:
    mat = sio.loadmat(str(priors_mat))
    if "Priors" not in mat:
        raise ValueError("{} missing variable 'Priors'".format(priors_mat))
    fc = mat_struct_get(mat["Priors"], "FC")
    return np.asarray(fc, dtype=np.float64)


def ridge_r2_raw_from_svd(maps: np.ndarray, priors_fc: np.ndarray, svd: Tuple[np.ndarray, np.ndarray, np.ndarray], lam: float) -> np.ndarray:
    u, s, vt = svd
    v = vt.T
    ut_y = u.T @ maps
    w = s / (s ** 2 + lam)
    beta = v @ (w[:, None] * ut_y)
    yhat = priors_fc @ beta
    sse = np.sum((maps - yhat) ** 2, axis=0)
    sst = np.sum((maps - np.mean(maps, axis=0, keepdims=True)) ** 2, axis=0)
    return 1.0 - sse / np.maximum(sst, np.finfo(np.float64).eps)


def ridge_r2_from_svd(maps: np.ndarray, priors_fc: np.ndarray, svd: Tuple[np.ndarray, np.ndarray, np.ndarray], lam: float) -> np.ndarray:
    r2 = ridge_r2_raw_from_svd(maps, priors_fc, svd, lam)
    # Published ICA-map NSI is ridge R2 against Priors.FC, explicitly clipped
    # into [0, 1]. Preserve raw negative values in diagnostics downstream.
    return np.where(np.isfinite(r2), np.clip(r2, 0.0, 1.0), 0.0)


def template_corr_details(maps: np.ndarray, priors_fc: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = maps - np.nanmean(maps, axis=0, keepdims=True)
    t = priors_fc - np.nanmean(priors_fc, axis=0, keepdims=True)
    num = t.T @ x
    den = np.sqrt(np.sum(t * t, axis=0, keepdims=True).T * np.sum(x * x, axis=0, keepdims=True))
    corr = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
    all_nan = np.all(~np.isfinite(corr), axis=0)
    safe_corr = np.where(np.isfinite(corr), corr, 0.0)
    best_idx = np.argmax(np.abs(safe_corr), axis=0)
    signed = safe_corr[best_idx, np.arange(safe_corr.shape[1])]
    max_abs = np.abs(signed)
    signed[all_nan] = np.nan
    max_abs[all_nan] = np.nan
    best_idx = best_idx.astype(float)
    best_idx[all_nan] = np.nan
    return max_abs, signed, best_idx + 1


def template_abs_corr(maps: np.ndarray, priors_fc: np.ndarray) -> np.ndarray:
    max_abs, _, _ = template_corr_details(maps, priors_fc)
    return max_abs


def load_neighbor_table(path: Path) -> np.ndarray:
    z = np.load(str(path), allow_pickle=False)
    if "neighbors" in z:
        return np.asarray(z["neighbors"])
    if len(z.files) == 1:
        return np.asarray(z[z.files[0]])
    raise ValueError("Could not find neighbors array in {}".format(path))


def build_cortex_adjacency(n_vertices: int, neighbor_path: Path) -> sp.csr_matrix:
    neighbors = load_neighbor_table(neighbor_path)
    if neighbors.shape[0] != n_vertices:
        raise ValueError(
            "Neighbor rows ({}) do not match cortical vertices ({})".format(neighbors.shape[0], n_vertices)
        )
    rows: List[int] = []
    cols: List[int] = []
    for i in range(n_vertices):
        for nbr in neighbors[i, 1:]:
            if np.isnan(nbr):
                continue
            j = int(nbr) - 1
            if 0 <= j < n_vertices:
                rows.append(i)
                cols.append(j)
    data = np.ones(len(rows), dtype=np.float64)
    w = sp.csr_matrix((data, (rows, cols)), shape=(n_vertices, n_vertices))
    return w.maximum(w.T).tocsr()


def morans_i(maps: np.ndarray, w: sp.csr_matrix) -> np.ndarray:
    n_vertices = maps.shape[0]
    centered = maps - np.mean(maps, axis=0, keepdims=True)
    weighted = w.dot(centered)
    num = np.sum(centered * weighted, axis=0)
    den = np.sum(centered * centered, axis=0)
    s0 = float(w.sum())
    return (n_vertices / s0) * (num / np.maximum(den, np.finfo(np.float64).eps))


def spectral_basis(w: sp.csr_matrix, kmax: int = 300) -> Tuple[np.ndarray, np.ndarray]:
    from scipy.sparse import csgraph
    from scipy.sparse import linalg as spla

    w = w.maximum(w.T).tocsr()
    degree = np.array(w.sum(axis=1)).ravel()
    degree[degree <= 0] = np.finfo(np.float64).eps
    dmh = sp.diags(1.0 / np.sqrt(degree))
    lap = sp.eye(w.shape[0], format="csr") - dmh @ w @ dmh
    n_conn = int(csgraph.connected_components(w, directed=False, return_labels=False))
    kreq = min(w.shape[0] - 1, int(kmax) + n_conn)
    evals = None
    evecs = None
    for sigma in (1e-6, 1e-5, 1e-4, 1e-3):
        try:
            evals, evecs = spla.eigsh(lap, k=kreq, sigma=sigma, which="LM")
            break
        except Exception:
            continue
    if evals is None:
        evals, evecs = spla.eigsh(lap + 1e-6 * sp.eye(w.shape[0], format="csr"), k=kreq, which="SA")
    order = np.argsort(evals)
    evals = evals[order][n_conn:]
    evecs = evecs[:, order][:, n_conn:]
    return evals, evecs


def spectral_slope_from_basis(maps: np.ndarray, basis: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    evals, evecs = basis
    x = maps - np.mean(maps, axis=0, keepdims=True)
    x = x / np.maximum(np.std(x, axis=0, ddof=1, keepdims=True), np.finfo(np.float64).eps)
    power = (evecs.T @ x) ** 2
    n_modes = evecs.shape[1]
    start = 5
    stop = n_modes - max(1, int(math.floor(0.10 * n_modes)))
    idx = np.arange(start, stop)
    log_freq = np.log(evals[idx] + 1e-12)
    slopes = np.full(maps.shape[1], np.nan, dtype=np.float64)
    design = np.column_stack([np.ones(idx.size), log_freq])
    for comp in range(maps.shape[1]):
        log_power = np.log(power[idx, comp] + 1e-12)
        beta, *_ = np.linalg.lstsq(design, log_power, rcond=None)
        slopes[comp] = beta[1]
    return slopes


def load_component_table(ted_dir: Path) -> pd.DataFrame:
    for metrics_tsv in (
        ted_dir / "ica_metrics.tsv",
        ted_dir / "desc-tedana_metrics.tsv",
        ted_dir / "comp_table.tsv",
    ):
        if not metrics_tsv.is_file():
            continue
        df = pd.read_csv(metrics_tsv, sep="\t")
        df = df.rename(
            columns={
                "component": "Component",
                "Kappa": "kappa",
                "Rho": "rho",
                "Classification": "classification",
            }
        )
        for col in ("Component", "kappa", "rho", "classification"):
            if col not in df.columns:
                raise ValueError("Missing required column '{}' in {}".format(col, metrics_tsv))
        df["component_id"] = df["Component"].map(component_id).astype(int)
        return df.sort_values("component_id").reset_index(drop=True)

    for decomp_json in (ted_dir / "ica_decomposition.json", ted_dir / "desc-ICA_decomposition.json"):
        if not decomp_json.is_file():
            continue
        obj = json.loads(decomp_json.read_text())
        rows: List[Dict[str, object]] = []
        for key, value in obj.items():
            if not str(key).startswith("ica_") or not isinstance(value, dict):
                continue
            row: Dict[str, object] = {"Component": str(key), "component_id": component_id(str(key))}
            for k, v in value.items():
                if k == "classification":
                    row[k] = str(v)
                elif np.isscalar(v) and not isinstance(v, (str, bytes)):
                    row[k] = v
            rows.append(row)
        if rows:
            df = pd.DataFrame(rows)
            for col in ("kappa", "rho", "classification"):
                if col not in df.columns:
                    raise ValueError("Missing required key '{}' in {}".format(col, decomp_json))
            return df.sort_values("component_id").reset_index(drop=True)

    raise FileNotFoundError("No component metrics table found in {}".format(ted_dir))


def run_parts(ted_dir: Path, root: Path) -> Dict[str, str]:
    rel = ted_dir.relative_to(root)
    parts = list(rel.parts)
    out = {
        "subject": parts[0] if len(parts) > 0 else "",
        "session": "",
        "run": "",
        "tedana_dir": str(ted_dir),
    }
    for part in parts:
        if part.startswith("session_"):
            out["session"] = part
        elif part.startswith("run_"):
            out["run"] = part
    return out


def read_component_maps(betas_cifti: Path, n_components: int, n_cortical_vertices: int) -> np.ndarray:
    img = nib.load(str(betas_cifti))
    data = np.asarray(img.dataobj, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError("Expected 2D CIFTI data, got shape {}".format(data.shape))
    if data.shape[1] == n_components:
        maps = data
    elif data.shape[0] == n_components:
        maps = data.T
    else:
        raise ValueError("CIFTI shape {} does not match {} components".format(data.shape, n_components))
    if maps.shape[0] < n_cortical_vertices:
        raise ValueError("CIFTI has fewer rows than Priors.FC")
    return maps[:n_cortical_vertices, :]


def iter_tedana_dirs(root: Path, subject_regex: Optional[str]) -> Iterable[Path]:
    pattern = re.compile(subject_regex) if subject_regex else None
    for betas in sorted(root.glob("*/func/rest/session_*/run_*/Tedana/betas_OC.dtseries.nii")):
        subject = betas.relative_to(root).parts[0]
        if pattern and not pattern.search(subject):
            continue
        ted_dir = betas.parent
        if (ted_dir / "ica_decomposition.json").is_file() or (ted_dir / "ica_metrics.tsv").is_file():
            yield ted_dir


def limit_runs_per_subject(ted_dirs: Sequence[Path], root: Path, max_runs: int) -> List[Path]:
    if max_runs <= 0:
        return list(ted_dirs)
    counts: Dict[str, int] = {}
    selected: List[Path] = []
    for ted_dir in ted_dirs:
        subject = ted_dir.relative_to(root).parts[0]
        count = counts.get(subject, 0)
        if count >= max_runs:
            continue
        selected.append(ted_dir)
        counts[subject] = count + 1
    return selected


def shard_path_for_run(ted_dir: Path, root: Path, shard_dir: Path) -> Path:
    rel = ted_dir.relative_to(root)
    stem = "__".join(rel.parts[:-1]) if rel.parts and rel.parts[-1] == "Tedana" else "__".join(rel.parts)
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)
    return shard_dir / "{}.csv".format(stem)


def extract_one_run(
    ted_dir: Path,
    root: Path,
    priors_fc: np.ndarray,
    svd: Tuple[np.ndarray, np.ndarray, np.ndarray],
    ridge_lambda: float,
    w: Optional[sp.csr_matrix],
    spectral: Optional[Tuple[np.ndarray, np.ndarray]],
) -> pd.DataFrame:
    df = load_component_table(ted_dir)
    maps = read_component_maps(ted_dir / "betas_OC.dtseries.nii", len(df), priors_fc.shape[0])
    run_info = run_parts(ted_dir, root)
    ridge_raw = ridge_r2_raw_from_svd(maps, priors_fc, svd, lam=ridge_lambda)
    nsi = np.where(np.isfinite(ridge_raw), np.clip(ridge_raw, 0.0, 1.0), 0.0)
    corr, signed_corr, best_prior = template_corr_details(maps, priors_fc)
    out = df.copy()
    for key, value in run_info.items():
        out[key] = value
    out["nsi"] = nsi
    out["nsi_fc_maxabs_corr"] = corr
    out["nsi_fc_best_signed_corr"] = signed_corr
    out["nsi_fc_best_prior"] = best_prior
    out["nsi_fc_ridge_r2_raw"] = ridge_raw
    out["nsi_fc_ridge_r2_pos"] = nsi
    out["nsi_ridge_negative_r2"] = np.isfinite(ridge_raw) & (ridge_raw < 0)
    out["nsi_ridge_nonfinite"] = ~np.isfinite(ridge_raw)
    out["template_abs_corr"] = corr
    out["map_mean"] = np.mean(maps, axis=0)
    out["map_std"] = np.std(maps, axis=0, ddof=1)
    out["map_abs_mean"] = np.mean(np.abs(maps), axis=0)
    out["map_abs_p95"] = np.percentile(np.abs(maps), 95, axis=0)
    if w is not None:
        out["morans_i"] = morans_i(maps, w)
    if spectral is not None:
        out["spectral_slope"] = spectral_slope_from_basis(maps, spectral)
    out["kappa"] = pd.to_numeric(out["kappa"], errors="coerce")
    out["rho"] = pd.to_numeric(out["rho"], errors="coerce")
    out["log10_kappa"] = np.log10(np.maximum(out["kappa"].to_numpy(dtype=float), np.finfo(float).eps))
    out["log10_rho"] = np.log10(np.maximum(out["rho"].to_numpy(dtype=float), np.finfo(float).eps))
    out["log10_kappa_rho"] = np.log10(
        out["kappa"].to_numpy(dtype=float) / np.maximum(out["rho"].to_numpy(dtype=float), np.finfo(float).eps)
    )
    return out


def extract_features(args: argparse.Namespace) -> pd.DataFrame:
    priors_fc = load_priors_fc(args.priors_mat)
    svd = np.linalg.svd(priors_fc, full_matrices=False)
    w = None
    spectral = None
    if args.compute_morans or args.compute_spectral_slope:
        print("[setup] building cortical adjacency", flush=True)
        w = build_cortex_adjacency(priors_fc.shape[0], args.neighbor_table)
    if args.compute_spectral_slope:
        print("[setup] computing spectral basis kmax={}".format(args.spectral_kmax), flush=True)
        spectral = spectral_basis(w, kmax=args.spectral_kmax)

    rows: List[pd.DataFrame] = []
    ted_dirs = list(iter_tedana_dirs(args.tms_root, args.subject_regex))
    ted_dirs = limit_runs_per_subject(ted_dirs, args.tms_root, args.max_runs_per_subject)
    if args.limit_runs:
        ted_dirs = ted_dirs[: args.limit_runs]
    args.shard_dir.mkdir(parents=True, exist_ok=True)
    for idx, ted_dir in enumerate(ted_dirs, start=1):
        shard = shard_path_for_run(ted_dir, args.tms_root, args.shard_dir)
        try:
            recompute = args.force_shards or not shard.is_file()
            if not recompute:
                out = pd.read_csv(shard)
                if args.compute_morans and "morans_i" not in out.columns:
                    recompute = True
                if args.compute_spectral_slope and "spectral_slope" not in out.columns:
                    recompute = True
            if recompute:
                out = extract_one_run(
                    ted_dir,
                    args.tms_root,
                    priors_fc,
                    svd,
                    args.ridge_lambda,
                    w if args.compute_morans else None,
                    spectral,
                )
                out.to_csv(shard, index=False)
            rows.append(out)
        except Exception as exc:
            print("[warn] skipping {}: {}".format(ted_dir, exc), flush=True)
        if idx % 25 == 0 or idx == len(ted_dirs):
            print("[extract] {}/{} runs".format(idx, len(ted_dirs)), flush=True)

    if not rows:
        raise RuntimeError("No rows extracted")
    all_rows = pd.concat(rows, ignore_index=True, sort=False)
    args.features_csv.parent.mkdir(parents=True, exist_ok=True)
    all_rows.to_csv(args.features_csv, index=False)
    return all_rows


def finite_xy(df: pd.DataFrame, features: Sequence[str], target: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    cols = list(features) + [target, "subject"]
    d = df.loc[:, cols].copy()
    d[target] = pd.to_numeric(d[target], errors="coerce")
    for col in features:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    keep = np.isfinite(d[target].to_numpy(dtype=float))
    d = d.loc[keep].reset_index(drop=True)
    return d.loc[:, features].to_numpy(dtype=float), d[target].to_numpy(dtype=float), d["subject"].astype(str).to_numpy(), d


def pearsonr_safe(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if y_true.size < 2 or np.std(y_true) <= 0 or np.std(y_pred) <= 0:
        return float("nan")
    return float(np.corrcoef(y_true, y_pred)[0, 1])


def model_specs(random_state: int, args: argparse.Namespace):
    specs = {
        "ridge_poly2": make_pipeline(
            SimpleImputer(strategy="median"),
            PolynomialFeatures(degree=2, include_bias=False),
            StandardScaler(),
            Ridge(alpha=1.0),
        ),
        "random_forest": make_pipeline(
            SimpleImputer(strategy="median"),
            RandomForestRegressor(
                n_estimators=args.rf_trees,
                min_samples_leaf=5,
                n_jobs=-1,
                random_state=random_state,
            ),
        ),
        "gradient_boosting": make_pipeline(
            SimpleImputer(strategy="median"),
            GradientBoostingRegressor(random_state=random_state),
        ),
    }
    wanted = [m.strip() for m in args.models.split(",") if m.strip()]
    unknown = sorted(set(wanted) - set(specs))
    if unknown:
        raise ValueError("Unknown model(s): {}".format(",".join(unknown)))
    return {name: specs[name] for name in wanted}


def evaluate_models(df: pd.DataFrame, args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame]:
    base_features = ["nsi"]
    feature_sets: Dict[str, List[str]] = {
        "nsi_only": base_features,
        "nsi_template": ["nsi", "template_abs_corr"],
        "nsi_map": ["nsi", "template_abs_corr", "map_std", "map_abs_mean", "map_abs_p95"],
    }
    if "morans_i" in df.columns:
        feature_sets["nsi_morans"] = ["nsi", "morans_i"]
        feature_sets["nsi_template_morans"] = ["nsi", "template_abs_corr", "morans_i"]
        feature_sets["nsi_map_morans"] = feature_sets["nsi_map"] + ["morans_i"]
    if "spectral_slope" in df.columns:
        feature_sets["nsi_slope"] = ["nsi", "spectral_slope"]
        feature_sets["nsi_template_slope"] = ["nsi", "template_abs_corr", "spectral_slope"]
        feature_sets["nsi_map_slope"] = feature_sets["nsi_map"] + ["spectral_slope"]
    if "morans_i" in df.columns and "spectral_slope" in df.columns:
        feature_sets["nsi_morans_slope"] = ["nsi", "morans_i", "spectral_slope"]
        feature_sets["nsi_template_morans_slope"] = [
            "nsi",
            "template_abs_corr",
            "morans_i",
            "spectral_slope",
        ]
        feature_sets["nsi_map_morans_slope"] = feature_sets["nsi_map"] + ["morans_i", "spectral_slope"]

    targets = ["log10_kappa_rho", "log10_kappa", "log10_rho"]
    metrics_rows: List[Dict[str, object]] = []
    pred_rows: List[pd.DataFrame] = []

    for target in targets:
        for feature_set_name, features in feature_sets.items():
            x, y, groups, usable = finite_xy(df, features, target)
            if len(np.unique(groups)) < 2:
                continue
            splitter = GroupShuffleSplit(
                n_splits=args.n_splits,
                train_size=args.train_fraction,
                random_state=args.random_seed,
            )
            for split_id, (train_idx, test_idx) in enumerate(splitter.split(x, y, groups), start=1):
                for model_name, model in model_specs(args.random_seed + split_id, args).items():
                    model.fit(x[train_idx], y[train_idx])
                    pred = np.asarray(model.predict(x[test_idx]), dtype=float)
                    truth = y[test_idx]
                    metrics_rows.append(
                        {
                            "target": target,
                            "feature_set": feature_set_name,
                            "model": model_name,
                            "split": split_id,
                            "n_train": int(train_idx.size),
                            "n_test": int(test_idx.size),
                            "n_train_subjects": int(np.unique(groups[train_idx]).size),
                            "n_test_subjects": int(np.unique(groups[test_idx]).size),
                            "r2": float(r2_score(truth, pred)),
                            "rmse": float(math.sqrt(mean_squared_error(truth, pred))),
                            "mae": float(mean_absolute_error(truth, pred)),
                            "pearson_r": pearsonr_safe(truth, pred),
                        }
                    )
                    if args.save_predictions:
                        p = usable.iloc[test_idx].copy()
                        p["target"] = target
                        p["feature_set"] = feature_set_name
                        p["model"] = model_name
                        p["split"] = split_id
                        p["y_true"] = truth
                        p["y_pred"] = pred
                        pred_rows.append(p)

    metrics = pd.DataFrame(metrics_rows)
    predictions = pd.concat(pred_rows, ignore_index=True, sort=False) if pred_rows else pd.DataFrame()
    args.metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.metrics_csv, index=False)
    if args.save_predictions:
        predictions.to_csv(args.predictions_csv, index=False)
    return metrics, predictions


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return metrics
    grouped = metrics.groupby(["target", "feature_set", "model"], as_index=False)
    summary = grouped.agg(
        r2_median=("r2", "median"),
        r2_mean=("r2", "mean"),
        rmse_median=("rmse", "median"),
        mae_median=("mae", "median"),
        pearson_r_median=("pearson_r", "median"),
        n_splits=("split", "count"),
    )
    return summary.sort_values(["target", "r2_median", "pearson_r_median"], ascending=[True, False, False])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tms-root", type=Path, required=True, help="Study data root for this optional training pilot")
    parser.add_argument("--priors-mat", type=Path, default=DEFAULT_PRIORS)
    parser.add_argument("--neighbor-table", type=Path, default=DEFAULT_NEIGHBORS)
    parser.add_argument("--features-csv", type=Path, default=Path("outputs/tms_pfm_component_features.csv"))
    parser.add_argument("--metrics-csv", type=Path, default=Path("outputs/tms_pfm_model_metrics.csv"))
    parser.add_argument("--predictions-csv", type=Path, default=Path("outputs/tms_pfm_model_predictions.csv"))
    parser.add_argument("--shard-dir", type=Path, default=Path("outputs/tms_pfm_feature_shards"))
    parser.add_argument("--force-shards", action="store_true")
    parser.add_argument("--reuse-features", action="store_true")
    parser.add_argument("--limit-runs", type=int, default=0)
    parser.add_argument("--max-runs-per-subject", type=int, default=0)
    parser.add_argument("--subject-regex", default="")
    parser.add_argument("--ridge-lambda", type=float, default=10.0)
    parser.add_argument("--compute-morans", action="store_true", default=True)
    parser.add_argument("--no-morans", dest="compute_morans", action="store_false")
    parser.add_argument("--compute-spectral-slope", action="store_true")
    parser.add_argument("--spectral-kmax", type=int, default=300)
    parser.add_argument("--n-splits", type=int, default=20)
    parser.add_argument("--train-fraction", type=float, default=0.5)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--models", default="ridge_poly2,random_forest,gradient_boosting")
    parser.add_argument("--rf-trees", type=int, default=100)
    parser.add_argument("--save-predictions", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.reuse_features and args.features_csv.is_file():
        df = pd.read_csv(args.features_csv)
    else:
        df = extract_features(args)
    metrics, _ = evaluate_models(df, args)
    summary = summarize(metrics)
    summary_path = args.metrics_csv.with_name(args.metrics_csv.stem + "_summary.csv")
    summary.to_csv(summary_path, index=False)
    print("\nExtracted rows: {} components from {} subjects".format(len(df), df["subject"].nunique()))
    print("Features: {}".format(args.features_csv))
    print("Metrics: {}".format(args.metrics_csv))
    print("Summary: {}".format(summary_path))
    if not summary.empty:
        print("\nTop model per target:")
        top = summary.sort_values(["target", "r2_median"], ascending=[True, False]).groupby("target").head(1)
        print(top.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate motion and grayplot QC reports for ME-fMRI processing stages."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mefmri_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np

from mefmri_motion import FD_METHOD, compute_fd_from_motion, filtered_motion_for_fd


REPORT_FILENAME = "motion_grayplots_qc.html"
GRAYPLOT_BASE_WIDTH_IN = 14.0
GRAYPLOT_FRAMES_PER_INCH = 180.0
GRAYPLOT_MAX_WIDTH_IN = 120.0


@dataclass
class RunInfo:
    session: int
    run: int
    run_dir: Path

    @property
    def tag(self) -> str:
        return f"S{self.session}_R{self.run}"

    @property
    def label(self) -> str:
        return f"S{self.session} R{self.run}"


@dataclass
class StageFile:
    tag: str
    path: Path


def parse_trailing_int(name: str, prefix: str) -> int | None:
    if not name.startswith(prefix):
        return None
    try:
        return int(name[len(prefix) :])
    except ValueError:
        return None


def discover_runs(subject_dir: Path, func_name: str, start_session: int) -> list[RunInfo]:
    root = subject_dir / "func" / func_name
    runs: list[RunInfo] = []
    for ses_dir in sorted(root.glob("session_*"), key=lambda p: parse_trailing_int(p.name, "session_") or 0):
        session = parse_trailing_int(ses_dir.name, "session_")
        if session is None or session < start_session:
            continue
        for run_dir in sorted(ses_dir.glob("run_*"), key=lambda p: parse_trailing_int(p.name, "run_") or 0):
            run = parse_trailing_int(run_dir.name, "run_")
            if run is None:
                continue
            runs.append(RunInfo(session, run, run_dir))
    return runs


def stage_rank(tag: str) -> tuple[int, str]:
    if tag == "OCME":
        return (0, tag)
    if tag.endswith("+MGTR"):
        return (90, tag)
    if "+MEICA" in tag:
        return (20, tag)
    if "+aCompCor" in tag:
        return (25, tag)
    if "+AROMA" in tag:
        return (30, tag)
    return (50, tag)


def discover_stage_files(run_dir: Path, func_prefix: str) -> list[StageFile]:
    stages: dict[str, Path] = {}
    for path in sorted(run_dir.glob(f"{func_prefix}_*.dtseries.nii")):
        name = path.name
        tag = name[len(func_prefix) + 1 : -len(".dtseries.nii")]
        if tag.endswith("_Betas") or tag.startswith("E"):
            continue
        if tag == "AVG_mcf":
            continue
        if tag == "OCME" or tag.startswith("OCME+"):
            stages[tag] = path
    canonical = [tag for tag in ("OCME", "OCME+MEICA", "OCME+MEICA+MGTR") if tag in stages]
    if len(canonical) >= 2:
        return [StageFile(tag, stages[tag]) for tag in canonical]
    return [StageFile(tag, stages[tag]) for tag in sorted(stages, key=stage_rank)]


def read_tr(run_dir: Path) -> float:
    tr_txt = run_dir / "TR.txt"
    if not tr_txt.exists():
        return 2.0
    try:
        return float(tr_txt.read_text().strip().split()[0])
    except Exception:
        return 2.0


def load_motion(run_dir: Path) -> np.ndarray | None:
    for name in ("MCF.par", "mcf.par"):
        path = run_dir / name
        if path.exists():
            motion = np.loadtxt(path, dtype=np.float64)
            if motion.ndim == 1:
                motion = motion.reshape(1, -1)
            if motion.shape[1] >= 6:
                return motion[:, :6]
    return None


def zscore_columns(data: np.ndarray) -> np.ndarray:
    data = np.asarray(data, dtype=np.float32)
    mu = np.nanmean(data, axis=0, keepdims=True)
    sig = np.nanstd(data, axis=0, keepdims=True)
    sig = np.where(sig > 1e-6, sig, 1.0)
    out = (data - mu) / sig
    out[~np.isfinite(out)] = 0
    return out.astype(np.float32, copy=False)


def zscore_vector(vec: np.ndarray) -> np.ndarray:
    vec = np.asarray(vec, dtype=np.float64)
    mu = np.nanmean(vec)
    sig = np.nanstd(vec)
    if not np.isfinite(sig) or sig <= 1e-8:
        return np.zeros(vec.shape, dtype=np.float64)
    out = (vec - mu) / sig
    out[~np.isfinite(out)] = 0
    return out


def dvars_from_z(z_time_by_gray: np.ndarray) -> np.ndarray:
    if z_time_by_gray.shape[0] == 0:
        return np.array([], dtype=np.float64)
    d = np.diff(z_time_by_gray, axis=0, prepend=z_time_by_gray[:1, :])
    return np.sqrt(np.mean(d * d, axis=1, dtype=np.float64))


def corr_1d(a: np.ndarray, b: np.ndarray) -> float:
    n = min(len(a), len(b))
    if n < 3:
        return float("nan")
    x = np.asarray(a[:n], dtype=np.float64)
    y = np.asarray(b[:n], dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(ok) < 3:
        return float("nan")
    x = x[ok] - np.mean(x[ok])
    y = y[ok] - np.mean(y[ok])
    denom = np.sqrt(np.sum(x * x) * np.sum(y * y))
    if denom <= 0:
        return float("nan")
    return float(np.sum(x * y) / denom)


def grayordinate_fd_corr(
    z_time_by_gray: np.ndarray,
    fd: np.ndarray | None,
    keep_mask: np.ndarray | None = None,
) -> tuple[float, float]:
    if fd is None or z_time_by_gray.shape[0] != len(fd):
        return float("nan"), float("nan")
    if keep_mask is not None:
        keep = np.asarray(keep_mask, dtype=bool)
        if keep.shape[0] != z_time_by_gray.shape[0] or np.count_nonzero(keep) < 3:
            return float("nan"), float("nan")
        z_time_by_gray = z_time_by_gray[keep, :]
        fd = fd[keep]
    fdz = zscore_vector(fd).astype(np.float32)
    z = z_time_by_gray - np.mean(z_time_by_gray, axis=0, keepdims=True)
    denom = np.sqrt(np.sum(z * z, axis=0) * np.sum(fdz * fdz))
    valid = denom > 1e-8
    if not np.any(valid):
        return float("nan"), float("nan")
    r = np.full(z.shape[1], np.nan, dtype=np.float64)
    r[valid] = np.sum(z[:, valid] * fdz[:, None], axis=0) / denom[valid]
    ar = np.abs(r[np.isfinite(r)])
    if ar.size == 0:
        return float("nan"), float("nan")
    return float(np.median(ar)), float(np.percentile(ar, 95))


def sample_grayordinates(img: nb.Cifti2Image, max_rows: int) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    axis = img.header.get_axis(1)
    structures = list(axis.iter_structures())
    total = int(img.shape[1])
    selected: list[np.ndarray] = []
    spans: list[tuple[str, int, int]] = []
    cursor = 0
    remaining = max_rows
    remaining_total = total

    for idx, (name, slc, _bm) in enumerate(structures):
        start = 0 if slc.start is None else int(slc.start)
        stop = total if slc.stop is None else int(slc.stop)
        size = max(0, stop - start)
        if size == 0:
            continue
        if idx == len(structures) - 1:
            n_take = remaining
        else:
            n_take = int(round(max_rows * size / max(total, 1)))
            if size >= 100:
                n_take = max(8, n_take)
            else:
                n_take = min(size, max(1, n_take))
            n_take = min(n_take, remaining)
        n_take = min(size, max(0, n_take))
        if n_take <= 0:
            remaining_total -= size
            continue
        if n_take >= size:
            idxs = np.arange(start, stop, dtype=np.int64)
        else:
            idxs = np.linspace(start, stop - 1, n_take, dtype=np.int64)
        selected.append(idxs)
        spans.append((short_structure_name(name), cursor, cursor + len(idxs)))
        cursor += len(idxs)
        remaining -= len(idxs)
        remaining_total -= size
        if remaining <= 0:
            break

    if not selected:
        idxs = np.linspace(0, total - 1, min(max_rows, total), dtype=np.int64)
        return idxs, [("grayordinates", 0, len(idxs))]
    return np.concatenate(selected), spans


def short_structure_name(name: str) -> str:
    return name.replace("CIFTI_STRUCTURE_", "").replace("_", " ").title()


def display_stage_label(tag: str) -> str:
    if tag == "OCME":
        return "OCME"
    if tag.endswith("+MGTR"):
        return "MGTR"
    if "+MEICA" in tag:
        return "MEICA"
    if "+aCompCor" in tag:
        return "aCompCor"
    if "+AROMA" in tag:
        return "AROMA"
    return tag.replace("OCME+", "")


def display_delta_label(a: str, b: str) -> str:
    return f"Delta {display_stage_label(b)} ({display_stage_label(b)} - {display_stage_label(a)})"


def stage_color(tag: str, idx: int = 0) -> str:
    if tag == "OCME":
        return "#4b5563"
    if "+MEICA" in tag:
        return "#2563eb"
    if tag.endswith("+MGTR"):
        return "#15803d"
    if "+aCompCor" in tag:
        return "#7c3aed"
    if "+AROMA" in tag:
        return "#c2410c"
    return plt.get_cmap("tab10").colors[idx % 10]


def load_sampled_dtseries(path: Path, sample_idx: np.ndarray) -> np.ndarray:
    img = nb.load(str(path))
    try:
        data = np.asanyarray(img.dataobj[:, sample_idx], dtype=np.float32)
    except Exception:
        full = np.asanyarray(img.dataobj, dtype=np.float32)
        data = full[:, sample_idx]
    data[~np.isfinite(data)] = 0
    return data


def plot_motion_grayplot_panel(
    out_png: Path,
    title: str,
    stages: list[str],
    stage_z: dict[str, np.ndarray],
    fd_raw: np.ndarray | None,
    fd_filtered: np.ndarray | None,
    motion: np.ndarray | None,
    dvars_by_stage: dict[str, np.ndarray],
    fd_threshold: float = 0.3,
    boundaries: list[int] | None = None,
) -> None:
    deltas = [(a, b) for a, b in zip(stages[:-1], stages[1:]) if a in stage_z and b in stage_z]
    n_rows = 3 + len(stages) + len(deltas)
    heights = [1.0, 0.5, 0.7] + [2.2 for _ in stages] + [1.6 for _ in deltas]
    n_frames = int(next(iter(stage_z.values())).shape[0])
    width_in = min(GRAYPLOT_MAX_WIDTH_IN, max(GRAYPLOT_BASE_WIDTH_IN, n_frames / GRAYPLOT_FRAMES_PER_INCH))
    fig, axes = plt.subplots(
        n_rows,
        1,
        figsize=(width_in, max(7, 1.2 + 1.25 * n_rows)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    x = np.arange(n_frames)

    ax = axes[0]
    if motion is not None:
        axis_colors = {"x": "#b91c1c", "y": "#1d4ed8", "z": "#15803d"}
        n_motion = min(n_frames, motion.shape[0])
        for idx, axis_name in enumerate(("x", "y", "z")):
            if motion.shape[1] > idx + 3:
                ax.plot(
                    x[:n_motion],
                    motion[:n_motion, idx + 3],
                    lw=0.75,
                    color=axis_colors[axis_name],
                    alpha=0.9,
                )
            if motion.shape[1] > idx:
                ax.plot(
                    x[:n_motion],
                    np.rad2deg(motion[:n_motion, idx]),
                    lw=0.75,
                    ls="--",
                    color=axis_colors[axis_name],
                    alpha=0.9,
                )
    ax.set_ylabel("motion")
    ax.set_title(title)
    ax.text(
        0.995,
        0.92,
        "translation (mm) / rotation (deg)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        color="#30343b",
    )
    ax.grid(True, axis="y", alpha=0.2)

    ax = axes[1]
    if fd_filtered is not None:
        ax.plot(x[: len(fd_filtered)], fd_filtered, color="#8b0000", lw=0.95)
        ax.axhline(fd_threshold, color="#8b0000", lw=0.7, ls="--", alpha=0.85)
        ax.text(
            0.01,
            fd_threshold,
            f"{fd_threshold:g} mm",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="bottom",
            fontsize=7,
            color="#8b0000",
        )
        censored = np.flatnonzero(np.asarray(fd_filtered) >= fd_threshold)
        if censored.size:
            ax.vlines(censored, ymin=0.0, ymax=fd_threshold, color="#8b0000", lw=0.35, alpha=0.25)
    ax.set_ylabel("FD")
    ax.grid(True, axis="y", alpha=0.2)

    ax = axes[2]
    dvars_items = []
    for idx, stage in enumerate(stages):
        dv = dvars_by_stage.get(stage)
        if dv is not None:
            dvars_items.append((idx, stage, dv))
    if dvars_items:
        offsets = np.linspace(3.2, -3.2, len(dvars_items)) if len(dvars_items) > 1 else np.array([0.0])
        for offset, (idx, stage, dv) in zip(offsets, dvars_items):
            color = stage_color(stage, idx)
            dv_z = np.clip(zscore_vector(dv), -2.8, 2.8) + offset
            ax.plot(x[: len(dv_z)], dv_z, lw=0.75, alpha=0.85, color=color)
            ax.text(
                1.002,
                offset,
                display_stage_label(stage),
                transform=ax.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=7,
                color=color,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.4},
                clip_on=False,
            )
        ax.set_yticks(offsets)
        ax.set_yticklabels([])
        ax.set_ylim(float(np.min(offsets) - 3.1), float(np.max(offsets) + 3.1))
    ax.set_ylabel("DVARS z\n+ offset")
    ax.grid(True, axis="y", alpha=0.2)

    stage_images = []
    delta_images = []

    row = 3
    for stage in stages:
        ax = axes[row]
        mat = np.clip(stage_z[stage].T, -2.5, 2.5)
        im = ax.imshow(mat, aspect="auto", interpolation="nearest", cmap="gray", vmin=-2.5, vmax=2.5)
        stage_images.append(im)
        ax.set_ylabel(display_stage_label(stage), fontsize=9, color=stage_color(stage, row - 3))
        ax.set_yticks([])
        row += 1

    for a, b in deltas:
        ax = axes[row]
        mat = np.clip((stage_z[b] - stage_z[a]).T, -2.5, 2.5)
        im = ax.imshow(mat, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        delta_images.append(im)
        delta_label = f"Delta {display_stage_label(b)}\n({display_stage_label(b)} - {display_stage_label(a)})"
        ax.set_ylabel("")
        ax.text(
            -0.012,
            0.5,
            delta_label,
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=8,
            clip_on=False,
        )
        ax.set_yticks([])
        row += 1

    if stage_images:
        cbar = fig.colorbar(stage_images[0], ax=list(axes[3 : 3 + len(stages)]), location="right", pad=0.006, fraction=0.016)
        cbar.set_label("z", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    if delta_images:
        cbar = fig.colorbar(delta_images[0], ax=list(axes[3 + len(stages) :]), location="right", pad=0.006, fraction=0.016)
        cbar.set_label("delta z", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    if boundaries:
        for ax in axes:
            for b in boundaries:
                ax.axvline(b - 0.5, color="#f0c36d", lw=0.7, alpha=0.9)

    axes[-1].set_xlabel("frame")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_metrics_tsv(rows: list[dict[str, object]], out_tsv: Path) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compute_qcfc_proxy(stage_run_samples: dict[str, list[np.ndarray]], run_fd_mean: list[float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    fd = np.asarray(run_fd_mean, dtype=np.float64)
    if len(fd) < 3 or np.nanstd(fd) <= 1e-8:
        return rows
    for stage, mats in stage_run_samples.items():
        if len(mats) != len(fd):
            continue
        edge_by_run = []
        for mat in mats:
            if mat.shape[0] < 5 or mat.shape[1] < 3:
                continue
            z = zscore_columns(mat)
            corr = np.corrcoef(z.T)
            iu = np.triu_indices(corr.shape[0], 1)
            edge_by_run.append(corr[iu].astype(np.float32))
        if len(edge_by_run) != len(fd):
            continue
        edges = np.vstack(edge_by_run)
        edge_r = []
        for col in range(edges.shape[1]):
            edge_r.append(corr_1d(edges[:, col], fd))
        edge_r = np.asarray(edge_r, dtype=np.float64)
        ar = np.abs(edge_r[np.isfinite(edge_r)])
        if ar.size == 0:
            continue
        rows.append(
            {
                "stage": stage,
                "n_runs": len(fd),
                "n_edges": int(ar.size),
                "median_abs_qcfc_r": float(np.median(ar)),
                "p95_abs_qcfc_r": float(np.percentile(ar, 95)),
                "frac_edges_abs_r_gt_0p2": float(np.mean(ar > 0.2)),
            }
        )
    return rows


def finite_float(value: object) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if np.isfinite(out) else float("nan")


def plot_motion_coupling_summary(out_png: Path, rows: list[dict[str, object]], fd_threshold: float) -> None:
    if not rows:
        return
    stages = sorted({str(row["stage"]) for row in rows}, key=stage_rank)
    run_keys = sorted({(int(row["session"]), int(row["run"])) for row in rows})
    run_labels = [f"S{s} R{r}" for s, r in run_keys]
    x = np.arange(len(run_keys), dtype=np.float64)
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(3, 1, figsize=(13.5, 8.2), sharex=True, constrained_layout=True)

    censor_by_run: list[float] = []
    max_fd_by_run: list[float] = []
    for run_key in run_keys:
        candidates = [row for row in rows if (int(row["session"]), int(row["run"])) == run_key]
        censor_by_run.append(finite_float(candidates[0].get("censor_fraction")) if candidates else float("nan"))
        max_fd_by_run.append(finite_float(candidates[0].get("max_fd")) if candidates else float("nan"))

    ax = axes[0]
    ax.bar(x, np.asarray(censor_by_run) * 100.0, color="#6f7c85", width=0.52, label=f"FD≥{fd_threshold:g}mm")
    ax2 = ax.twinx()
    ax2.plot(x, max_fd_by_run, color="#8b0000", marker="o", lw=1.3, label="max FD")
    ax.set_ylabel("censored %")
    ax2.set_ylabel("max FD mm")
    ax.set_title("Motion Censoring Burden")
    ax.grid(True, axis="y", alpha=0.2)
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, loc="upper right", fontsize=8)

    ax = axes[1]
    for idx, stage in enumerate(stages):
        vals_all = []
        vals_keep = []
        for run_key in run_keys:
            match = next(
                (
                    row
                    for row in rows
                    if str(row["stage"]) == stage and (int(row["session"]), int(row["run"])) == run_key
                ),
                None,
            )
            vals_all.append(finite_float(match.get("abs_corr_dvars_fd")) if match else float("nan"))
            vals_keep.append(finite_float(match.get("abs_corr_dvars_fd_retained")) if match else float("nan"))
        color = colors[idx % len(colors)]
        ax.plot(x, vals_all, color=color, lw=0.9, alpha=0.28, ls="--")
        ax.plot(x, vals_keep, color=color, marker="o", lw=1.6, label=display_stage_label(stage))
    ax.set_ylabel("abs r")
    ax.set_title("DVARS-FD Coupling After FD Censoring")
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(loc="upper right", ncol=min(4, len(stages)), fontsize=8)

    ax = axes[2]
    for idx, stage in enumerate(stages):
        vals_all = []
        vals_keep = []
        for run_key in run_keys:
            match = next(
                (
                    row
                    for row in rows
                    if str(row["stage"]) == stage and (int(row["session"]), int(row["run"])) == run_key
                ),
                None,
            )
            vals_all.append(finite_float(match.get("p95_abs_grayordinate_fd_r")) if match else float("nan"))
            vals_keep.append(finite_float(match.get("p95_abs_grayordinate_fd_r_retained")) if match else float("nan"))
        color = colors[idx % len(colors)]
        ax.plot(x, vals_all, color=color, lw=0.9, alpha=0.28, ls="--")
        ax.plot(x, vals_keep, color=color, marker="o", lw=1.6, label=display_stage_label(stage))
    ax.set_ylabel("p95 abs r")
    ax.set_title("Grayordinate-FD Coupling After FD Censoring")
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", alpha=0.2)
    ax.set_xticks(x)
    ax.set_xticklabels(run_labels)
    ax.set_xlabel("run")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_qcfc_proxy_summary(out_png: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    rows = sorted(rows, key=lambda row: stage_rank(str(row["stage"])))
    labels = [display_stage_label(str(row["stage"])) for row in rows]
    x = np.arange(len(rows), dtype=np.float64)
    median_abs = np.asarray([finite_float(row.get("median_abs_qcfc_r")) for row in rows], dtype=np.float64)
    p95_abs = np.asarray([finite_float(row.get("p95_abs_qcfc_r")) for row in rows], dtype=np.float64)
    frac = np.asarray([finite_float(row.get("frac_edges_abs_r_gt_0p2")) for row in rows], dtype=np.float64)

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 4.8), constrained_layout=True)
    colors = ["#6f7c85", "#2f80ed", "#2e7d32", "#8b5cf6", "#c2410c"]
    ax.bar(x, median_abs, color=[colors[i % len(colors)] for i in range(len(rows))], width=0.58, alpha=0.86, label="median abs QC-FC r")
    ax.scatter(x, p95_abs, color="#111827", s=42, marker="D", label="p95 abs QC-FC r", zorder=3)
    for xi, f in zip(x, frac):
        if np.isfinite(f):
            ax.text(xi, 0.02, f"{100.0 * f:.0f}%\n>|.2|", ha="center", va="bottom", fontsize=8, color="white")
    ax.set_ylabel("abs run-level QC-FC r")
    ax.set_title("Exploratory QC-FC Motion Coupling Proxy")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, max(1.0, float(np.nanmax(p95_abs)) * 1.08 if np.any(np.isfinite(p95_abs)) else 1.0))
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(loc="upper left", fontsize=8)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def generate_report(
    subject_dir: Path,
    func_name: str,
    func_prefix: str,
    start_session: int,
    out_dir: Path,
    max_grayplot_rows: int,
    qcfc_sample_size: int,
    fd_threshold: float,
) -> None:
    runs = discover_runs(subject_dir, func_name, start_session)
    if not runs:
        raise SystemExit(f"No runs found under {subject_dir / 'func' / func_name}")

    out_dir.mkdir(parents=True, exist_ok=True)
    reportlets = out_dir / "reportlets"
    reportlets.mkdir(parents=True, exist_ok=True)

    first_stage = None
    first_img = None
    run_stages: dict[str, list[StageFile]] = {}
    for run in runs:
        stages = discover_stage_files(run.run_dir, func_prefix)
        if stages:
            run_stages[run.tag] = stages
            if first_stage is None:
                first_stage = stages[0]
                first_img = nb.load(str(first_stage.path))
    if first_img is None:
        raise SystemExit(f"No {func_prefix}_*.dtseries.nii stage files found")

    sample_idx, structure_spans = sample_grayordinates(first_img, max_grayplot_rows)
    all_stage_tags = sorted({stage.tag for stages in run_stages.values() for stage in stages}, key=stage_rank)
    run_metrics: list[dict[str, object]] = []
    run_cards: list[str] = []
    run_links: list[str] = []
    concat_stage_blocks: dict[str, list[np.ndarray]] = {tag: [] for tag in all_stage_tags}
    concat_dvars: dict[str, list[np.ndarray]] = {tag: [] for tag in all_stage_tags}
    concat_fd_raw: list[np.ndarray] = []
    concat_fd: list[np.ndarray] = []
    concat_motion: list[np.ndarray] = []
    boundaries: list[int] = []
    run_fd_mean: list[float] = []

    cursor = 0
    for run in runs:
        stages = run_stages.get(run.tag, [])
        if not stages:
            continue
        tr = read_tr(run.run_dir)
        motion = load_motion(run.run_dir)
        fd_raw = compute_fd_from_motion(motion, tr)
        filtered_motion, fd_filter = filtered_motion_for_fd(motion, tr)
        fd = compute_fd_from_motion(filtered_motion, tr)
        if fd is not None:
            keep_mask = np.asarray(fd) < fd_threshold
            n_censored = int(np.count_nonzero(~keep_mask))
            censor_fraction = float(n_censored / max(1, len(fd)))
        else:
            keep_mask = None
            n_censored = 0
            censor_fraction = float("nan")
        stage_z: dict[str, np.ndarray] = {}
        dvars_by_stage: dict[str, np.ndarray] = {}

        for stage in stages:
            sampled = load_sampled_dtseries(stage.path, sample_idx)
            z = zscore_columns(sampled)
            stage_z[stage.tag] = z
            dvars_by_stage[stage.tag] = dvars_from_z(z)

            gsig = np.mean(z, axis=1)
            med_abs, p95_abs = grayordinate_fd_corr(z, fd)
            med_abs_keep, p95_abs_keep = grayordinate_fd_corr(z, fd, keep_mask)
            run_metrics.append(
                {
                    "session": run.session,
                    "run": run.run,
                    "stage": stage.tag,
                    "n_timepoints": z.shape[0],
                    "fd_threshold_mm": fd_threshold,
                    "fd_filter": fd_filter,
                    "n_fd_censored": n_censored,
                    "censor_fraction": censor_fraction,
                    "n_retained": int(np.count_nonzero(keep_mask)) if keep_mask is not None else "",
                    "mean_raw_fd": float(np.nanmean(fd_raw)) if fd_raw is not None else float("nan"),
                    "max_raw_fd": float(np.nanmax(fd_raw)) if fd_raw is not None else float("nan"),
                    "mean_fd": float(np.nanmean(fd)) if fd is not None else float("nan"),
                    "max_fd": float(np.nanmax(fd)) if fd is not None else float("nan"),
                    "mean_dvars": float(np.nanmean(dvars_by_stage[stage.tag])),
                    "corr_dvars_fd": corr_1d(dvars_by_stage[stage.tag], fd) if fd is not None else float("nan"),
                    "abs_corr_dvars_fd": abs(corr_1d(dvars_by_stage[stage.tag], fd)) if fd is not None else float("nan"),
                    "corr_dvars_fd_retained": corr_1d(dvars_by_stage[stage.tag][keep_mask], fd[keep_mask]) if fd is not None and keep_mask is not None and np.count_nonzero(keep_mask) >= 3 else float("nan"),
                    "abs_corr_dvars_fd_retained": abs(corr_1d(dvars_by_stage[stage.tag][keep_mask], fd[keep_mask])) if fd is not None and keep_mask is not None and np.count_nonzero(keep_mask) >= 3 else float("nan"),
                    "corr_global_signal_fd": corr_1d(gsig, fd) if fd is not None else float("nan"),
                    "corr_global_signal_fd_retained": corr_1d(gsig[keep_mask], fd[keep_mask]) if fd is not None and keep_mask is not None and np.count_nonzero(keep_mask) >= 3 else float("nan"),
                    "median_abs_grayordinate_fd_r": med_abs,
                    "p95_abs_grayordinate_fd_r": p95_abs,
                    "median_abs_grayordinate_fd_r_retained": med_abs_keep,
                    "p95_abs_grayordinate_fd_r_retained": p95_abs_keep,
                    "dtseries": str(stage.path),
                }
            )

        common_stages = [tag for tag in all_stage_tags if tag in stage_z]
        run_png = reportlets / f"subj_desc-motionGrayplot_{run.tag}.png"
        plot_motion_grayplot_panel(
            run_png,
            f"{subject_dir.name} {func_name} {run.label}",
            common_stages,
            stage_z,
            fd_raw,
            fd,
            motion,
            dvars_by_stage,
            fd_threshold=fd_threshold,
        )
        run_cards.append(
            "<section class='run'>"
            f"<h3>{html.escape(run.label)}</h3>"
            f"<a href='reportlets/{html.escape(run_png.name)}'>"
            f"<img src='reportlets/{html.escape(run_png.name)}' alt='Motion grayplot {html.escape(run.tag)}'>"
            "</a>"
            "</section>"
        )
        run_links.append(f"<a href='reportlets/{html.escape(run_png.name)}'>{html.escape(run.label)}</a>")

        n_tp = next(iter(stage_z.values())).shape[0]
        for tag in all_stage_tags:
            if tag in stage_z:
                concat_stage_blocks[tag].append(stage_z[tag])
                concat_dvars[tag].append(dvars_by_stage[tag])
        if fd_raw is not None:
            concat_fd_raw.append(fd_raw)
        else:
            concat_fd_raw.append(np.full(n_tp, np.nan))
        if fd is not None:
            concat_fd.append(fd)
            run_fd_mean.append(float(np.nanmean(fd)))
        else:
            concat_fd.append(np.full(n_tp, np.nan))
            run_fd_mean.append(float("nan"))
        if motion is not None:
            concat_motion.append(motion)
        else:
            concat_motion.append(np.full((n_tp, 6), np.nan))
        cursor += n_tp
        boundaries.append(cursor)

    common_concat_stages = [tag for tag in all_stage_tags if len(concat_stage_blocks[tag]) == len(run_fd_mean)]
    concat_stage_z = {tag: np.vstack(concat_stage_blocks[tag]) for tag in common_concat_stages}
    concat_dvars_z = {tag: np.concatenate(concat_dvars[tag]) for tag in common_concat_stages}
    concat_fd_raw_vec = np.concatenate(concat_fd_raw) if concat_fd_raw else None
    concat_fd_vec = np.concatenate(concat_fd) if concat_fd else None
    concat_motion_arr = np.vstack(concat_motion) if concat_motion else None
    subject_png = reportlets / "subj_desc-motionGrayplots_stacked.png"
    if concat_stage_z:
        plot_motion_grayplot_panel(
            subject_png,
            f"{subject_dir.name} {func_name} concatenated runs",
            common_concat_stages,
            concat_stage_z,
            concat_fd_raw_vec,
            concat_fd_vec,
            concat_motion_arr,
            concat_dvars_z,
            fd_threshold=fd_threshold,
            boundaries=boundaries[:-1],
        )

    write_metrics_tsv(run_metrics, out_dir / "motion_grayplot_metrics.tsv")

    coupling_png = reportlets / "subj_desc-motionCouplingSummary.png"
    plot_motion_coupling_summary(coupling_png, run_metrics, fd_threshold)

    subject_img = (
        f"<section><h2>Stacked Grayplots Across Runs</h2><div class='scroll-x'><img class='wide-img' style='width: {max(100.0, 100.0 * len(run_cards) / 4.0):.1f}%;' src='reportlets/{html.escape(subject_png.name)}' alt='Stacked motion grayplots'></div></section>"
        if subject_png.exists()
        else ""
    )
    coupling_img = (
        f"<section><h2>Motion Coupling Summary</h2><img src='reportlets/{html.escape(coupling_png.name)}' alt='Motion coupling summary'></section>"
        if coupling_png.exists()
        else ""
    )
    interpretation_note = (
        "<section><h2>Motion Coupling Interpretation</h2>"
        "<p><strong>abs r</strong> is the absolute correlation between a motion estimate and a signal-change metric. "
        "It ignores sign because both positive and negative coupling mean the data still vary systematically with motion. "
        "Values closer to zero are better.</p>"
        "<p><strong>p95 abs r</strong> is the 95th percentile of absolute grayordinate-FD correlations. "
        "It summarizes the high-coupling tail rather than the average grayordinate. Lower values indicate fewer strongly motion-coupled grayordinates.</p>"
        "<p>Dashed lines in the summary use all frames. Solid lines use retained frames after FD censoring. "
        "Motion is bandstop-filtered from 0.2 Hz to Nyquist minus 0.019 Hz only for runs with TR &lt; 2 s, "
        "using the same FD calculation as CIFTI concatenation.</p>"
        "</section>"
    )
    full_size_links = (
        "<section><h2>Full-Size Run Grayplots</h2><div class='link-row'>"
        + " ".join(run_links)
        + "</div></section>"
        if run_links
        else ""
    )
    report = out_dir / REPORT_FILENAME
    report.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(subject_dir.name)} {html.escape(func_name)} Motion Grayplots QC</title>
<style>
body {{ margin: 0; font-family: Arial, sans-serif; color: #1f2328; background: #f6f7f9; }}
header {{ padding: 22px 28px; background: #ffffff; border-bottom: 1px solid #d8dee4; }}
main {{ max-width: 1280px; margin: 0 auto; padding: 22px; }}
section {{ background: #ffffff; border: 1px solid #d8dee4; border-radius: 6px; padding: 16px; margin-bottom: 18px; }}
h1, h2, h3 {{ margin: 0 0 12px; }}
img {{ max-width: 100%; height: auto; display: block; margin: 10px 0; border: 1px solid #d8dee4; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th, td {{ border-bottom: 1px solid #d8dee4; padding: 7px 9px; text-align: left; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(520px, 1fr)); gap: 18px; }}
.scroll-x {{ overflow-x: auto; overflow-y: hidden; padding-bottom: 8px; }}
.scroll-x .wide-img {{ max-width: none; }}
.run a {{ display: block; color: inherit; text-decoration: none; }}
.link-row {{ display: flex; flex-wrap: wrap; gap: 10px; }}
.link-row a {{ color: #0969da; text-decoration: none; border: 1px solid #d8dee4; border-radius: 4px; padding: 6px 9px; background: #ffffff; }}
.link-row a:hover {{ text-decoration: underline; }}
</style>
</head>
<body>
<header>
<h1>{html.escape(subject_dir.name)} {html.escape(func_name)} Motion &amp; Grayplots QC</h1>
<div>Stacked grayplots, removed-signal deltas, motion, FD, DVARS, FD≥{fd_threshold:g} mm censor markers, and motion-coupling summaries.</div>
</header>
<main>
{subject_img}
<div class="grid">
{''.join(run_cards)}
</div>
{coupling_img}
{interpretation_note}
{full_size_links}
<section><h2>Detailed Metrics</h2><div><a href="motion_grayplot_metrics.tsv">motion_grayplot_metrics.tsv</a> contains run/stage motion coupling values.</div></section>
</main>
</body>
</html>
""",
        encoding="utf-8",
    )

    payload = {
        "subject_dir": str(subject_dir),
        "func_name": func_name,
        "func_prefix": func_prefix,
        "start_session": start_session,
        "stages": all_stage_tags,
        "n_runs": len(run_cards),
        "sampled_grayordinates": int(len(sample_idx)),
        "fd_threshold_mm": fd_threshold,
        "fd_method": FD_METHOD,
        "fd_keep_comparison": "<",
        "report": str(report),
    }
    (out_dir / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(report)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject-dir", required=True, type=Path)
    parser.add_argument("--func-name", default="rest")
    parser.add_argument("--func-prefix", default="Rest")
    parser.add_argument("--start-session", type=int, default=1)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--max-grayplot-rows", type=int, default=1200)
    parser.add_argument("--qcfc-sample-size", type=int, default=300)
    parser.add_argument("--fd-threshold", type=float, default=0.3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    subject_dir = args.subject_dir.resolve()
    qc_dirname = os.environ.get("QC_OUTPUT_DIRNAME", "QC")
    out_dir = (args.out_dir or subject_dir / "func" / args.func_name / qc_dirname / "MotionGrayPlots").resolve()
    generate_report(
        subject_dir,
        args.func_name,
        args.func_prefix,
        args.start_session,
        out_dir,
        args.max_grayplot_rows,
        args.qcfc_sample_size,
        args.fd_threshold,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

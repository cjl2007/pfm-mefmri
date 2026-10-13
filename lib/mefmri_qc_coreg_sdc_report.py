#!/usr/bin/env python3
"""Experimental fMRIPrep-like QC report for ME-fMRI SDC/coregistration.

This is intentionally standalone: it reads existing pipeline outputs and writes
static PNG/SVG/GIF/HTML reportlets without changing preprocessing outputs.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mefmri_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from PIL import Image


COREG_REPORT_FILENAME = "functional_coregistration_qc.html"
FIELDMAP_REPORT_FILENAME = "distortion_correction_qc.html"


@dataclass
class RunInfo:
    session: int
    run: int
    run_dir: Path
    method: str
    distorted: Path | None
    unwarped: Path | None
    corrected_acpc: Path | None
    displacement: Path | None
    pe_axis: str

    @property
    def tag(self) -> str:
        return f"S{self.session}_R{self.run}"

    @property
    def label(self) -> str:
        return f"session {self.session}, run {self.run}"


def _read_text(path: Path, default: str = "") -> str:
    try:
        return path.read_text().strip()
    except OSError:
        return default


def discover_runs(subject_dir: Path, func_name: str, xfms_name: str, distortion_mode: str = "auto") -> list[RunInfo]:
    runs: list[RunInfo] = []
    func_root = subject_dir / "func" / func_name
    for run_dir in sorted(func_root.glob("session_*/run_*")):
        match = re.search(r"session_(\d+)/run_(\d+)$", str(run_dir))
        if not match:
            continue
        session, run = (int(match.group(1)), int(match.group(2)))
        ref_dir = run_dir / "MEDIC" / "reference"
        distorted = ref_dir / "SBref_distorted.nii.gz"
        unwarped = ref_dir / "SBref_unwarped.nii.gz"
        use_medic = distortion_mode in {"auto", "medic"} and distorted.exists() and unwarped.exists()
        if use_medic:
            disp_txt = run_dir / "MEDICDisplacementMap.txt"
            disp_path = Path(_read_text(disp_txt)) if disp_txt.exists() else None
            if disp_path and not disp_path.exists():
                disp_path = None
            pe_axis = _read_text(run_dir / "MEDICPhaseEncodingAxis.txt", _read_text(run_dir / "PE.txt", "n/a"))
            runs.append(
                RunInfo(
                    session=session,
                    run=run,
                    run_dir=run_dir,
                    method="medic",
                    distorted=distorted,
                    unwarped=unwarped,
                    corrected_acpc=None,
                    displacement=disp_path,
                    pe_axis=pe_axis,
                )
            )
            continue

        if distortion_mode == "medic":
            continue

        sbref = run_dir / "SBref.nii.gz"
        if not sbref.exists():
            continue
        qc_dirname = os.environ.get("QC_OUTPUT_DIRNAME", "QC")
        coreg_dir = subject_dir / "func" / func_name / qc_dirname / "Coreg"
        if not coreg_dir.exists():
            coreg_dir = subject_dir / "func" / func_name / "qa" / "Coreg"
        scan_specific = coreg_dir / f"SBref2acpc_EpiReg+BBR_ScanSpecificFM_S{session}_R{run}.nii.gz"
        avgfm = coreg_dir / f"SBref2acpc_EpiReg+BBR_AvgFM_S{session}_R{run}.nii.gz"
        corrected_acpc = scan_specific if scan_specific.exists() else avgfm if avgfm.exists() else None
        if corrected_acpc is None:
            continue
        pe_axis = _read_text(run_dir / "PE.txt", "n/a")
        fmap = subject_dir / "func" / func_name / "field_maps" / "AllFMs" / f"FM_rads_acpc_S{session}_R{run}.nii.gz"
        method = "fieldmap" if distortion_mode == "auto" else distortion_mode
        runs.append(
            RunInfo(
                session=session,
                run=run,
                run_dir=run_dir,
                method=method,
                distorted=sbref,
                unwarped=None,
                corrected_acpc=corrected_acpc,
                displacement=fmap if fmap.exists() else None,
                pe_axis=pe_axis,
            )
        )
    return runs


def load_3d(path: Path, volume: int | None = None) -> np.ndarray:
    img = nb.load(str(path))
    data = np.asanyarray(img.dataobj)
    if data.ndim == 4:
        data = data[..., 0 if volume is None else volume]
    data = np.asarray(data, dtype=np.float32)
    data[~np.isfinite(data)] = 0
    return data


def resample_run_unwarped_sbref_to_acpc(run: RunInfo, reference: Path, work_dir: Path) -> Path | None:
    if run.method != "medic" and run.corrected_acpc is not None:
        return run.corrected_acpc
    if run.unwarped is None:
        return None
    sbref = Path(_read_text(run.run_dir / "MEDICUnwarpedSBRef.txt", str(run.unwarped)))
    target = Path(_read_text(run.run_dir / "IntermediateCoregTarget.txt"))
    warp = Path(_read_text(run.run_dir / "Intermediate2ACPCWarp.txt"))
    if not sbref.exists() or not target.exists() or not warp.exists():
        return None
    out_file = work_dir / f"subj_desc-coreg_{run.tag}_unwarpedSBRef_acpc.nii.gz"
    mat_file = work_dir / f"subj_desc-coreg_{run.tag}_unwarpedSBRef_toIntermediate.mat"
    target_file = work_dir / f"subj_desc-coreg_{run.tag}_unwarpedSBRef_toIntermediate.nii.gz"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "flirt",
            "-dof",
            "6",
            "-in",
            str(sbref),
            "-ref",
            str(target),
            "-omat",
            str(mat_file),
            "-out",
            str(target_file),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "applywarp",
            "--interp=spline",
            f"--in={sbref}",
            f"--premat={mat_file}",
            f"--ref={reference}",
            f"--out={out_file}",
            f"--warp={warp}",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return out_file


def resample_medic_sdc_pair_to_acpc(run: RunInfo, reference: Path, work_dir: Path) -> tuple[Path, Path] | None:
    if run.distorted is None or run.unwarped is None:
        return None
    unwarped = Path(_read_text(run.run_dir / "MEDICUnwarpedSBRef.txt", str(run.unwarped)))
    target = Path(_read_text(run.run_dir / "IntermediateCoregTarget.txt"))
    warp = Path(_read_text(run.run_dir / "Intermediate2ACPCWarp.txt"))
    if not run.distorted.exists() or not unwarped.exists() or not target.exists() or not warp.exists():
        return None

    work_dir.mkdir(parents=True, exist_ok=True)
    mat_file = work_dir / f"subj_desc-sdc_{run.tag}_unwarpedSBRef_toIntermediate.mat"
    target_file = work_dir / f"subj_desc-sdc_{run.tag}_unwarpedSBRef_toIntermediate.nii.gz"
    distorted_out = work_dir / f"subj_desc-sdc_{run.tag}_distorted_acpc.nii.gz"
    unwarped_out = work_dir / f"subj_desc-sdc_{run.tag}_unwarped_acpc.nii.gz"

    subprocess.run(
        [
            "flirt",
            "-dof",
            "6",
            "-in",
            str(unwarped),
            "-ref",
            str(target),
            "-omat",
            str(mat_file),
            "-out",
            str(target_file),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    for in_file, out_file in ((run.distorted, distorted_out), (unwarped, unwarped_out)):
        subprocess.run(
            [
                "applywarp",
                "--interp=spline",
                f"--in={in_file}",
                f"--premat={mat_file}",
                f"--ref={reference}",
                f"--out={out_file}",
                f"--warp={warp}",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    return distorted_out, unwarped_out


def resample_standard_sdc_pair_to_acpc(
    run: RunInfo,
    subject_dir: Path,
    func_name: str,
    xfms_name: str,
    reference: Path,
    work_dir: Path,
) -> tuple[Path, Path] | None:
    if run.distorted is None or run.corrected_acpc is None:
        return None
    xfm_dir = subject_dir / "func" / "xfms" / xfms_name
    qc_dirname = os.environ.get("QC_OUTPUT_DIRNAME", "QC")
    coreg_dir = subject_dir / "func" / func_name / qc_dirname / "Coreg"
    if not coreg_dir.exists():
        coreg_dir = subject_dir / "func" / func_name / "qa" / "Coreg"
    run_to_avg = coreg_dir / f"SBref2AvgSBref_S{run.session}_R{run.run}.mat"
    avg_to_acpc = xfm_dir / "AvgSBref2acpc_EpiReg+BBR_full.mat"
    if not run.distorted.exists() or not run.corrected_acpc.exists() or not run_to_avg.exists() or not avg_to_acpc.exists():
        return None

    work_dir.mkdir(parents=True, exist_ok=True)
    concat_mat = work_dir / f"subj_desc-sdc_{run.tag}_distorted_toAcpc_affine.mat"
    distorted_out = work_dir / f"subj_desc-sdc_{run.tag}_distorted_acpc_affine.nii.gz"
    subprocess.run(
        ["convert_xfm", "-omat", str(concat_mat), "-concat", str(avg_to_acpc), str(run_to_avg)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "flirt",
            "-dof",
            "6",
            "-interp",
            "spline",
            "-in",
            str(run.distorted),
            "-ref",
            str(reference),
            "-out",
            str(distorted_out),
            "-applyxfm",
            "-init",
            str(concat_mat),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return distorted_out, run.corrected_acpc


def robust_limits(data: np.ndarray, symmetric: bool = False) -> tuple[float, float]:
    values = data[np.isfinite(data)]
    if values.size == 0:
        return 0.0, 1.0
    if symmetric:
        lim = float(np.percentile(np.abs(values), 99))
        lim = max(lim, 1e-3)
        return -lim, lim
    nz = values[np.abs(values) > 1e-8]
    if nz.size == 0:
        nz = values
    lo, hi = np.percentile(nz, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.min(values)), float(np.max(values))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def center_slices(data: np.ndarray, mask: np.ndarray | None = None) -> tuple[int, int, int]:
    if mask is None:
        finite = np.isfinite(data)
        nz = finite & (np.abs(data) > np.percentile(np.abs(data[finite]), 20) if finite.any() else False)
        mask = nz
    coords = np.argwhere(mask)
    if coords.size == 0:
        return tuple(int(s // 2) for s in data.shape[:3])
    return tuple(int(round(v)) for v in coords.mean(axis=0))


def off_midline_sagittal_slices(
    data: np.ndarray,
    contour: np.ndarray | None,
    fraction: float = 0.07,
) -> tuple[int, int, int]:
    slices = list(center_slices(data, contour))
    if contour is None or contour.shape[0] < 3:
        return tuple(slices)

    offset = max(2, int(round(contour.shape[0] * fraction)))
    center = slices[0]
    candidates = [idx for idx in (center - offset, center + offset) if 0 <= idx < contour.shape[0]]
    if not candidates:
        return tuple(slices)

    slices[0] = max(candidates, key=lambda idx: int(np.count_nonzero(contour[idx, :, :])))
    return tuple(slices)


def slice2d(data: np.ndarray, axis: int, index: int) -> np.ndarray:
    if axis == 0:
        out = data[index, :, :]
    elif axis == 1:
        out = data[:, index, :]
    else:
        out = data[:, :, index]
    return np.flipud(np.rot90(out))


def make_mask(data: np.ndarray) -> np.ndarray:
    finite = np.isfinite(data)
    vals = data[finite]
    if vals.size == 0:
        return np.zeros(data.shape, dtype=bool)
    thr = np.percentile(vals, 35)
    return finite & (data > thr)


def plot_rows(
    rows: list[tuple[str, np.ndarray, str, bool]],
    out_png: Path,
    title: str,
    contour: np.ndarray | None = None,
    contour_color: str = "#d62728",
    slices: tuple[int, int, int] | None = None,
) -> None:
    if slices is None:
        slices = center_slices(rows[0][1], contour)
    axes_names = ["Sagittal", "Coronal", "Axial"]
    fig, axs = plt.subplots(
        len(rows),
        3,
        figsize=(10.5, 3.0 * len(rows)),
        facecolor="white",
        constrained_layout=True,
    )
    if len(rows) == 1:
        axs = np.asarray([axs])
    fig.suptitle(title, fontsize=14, fontweight="bold")

    for row_idx, (label, data, cmap, symmetric) in enumerate(rows):
        vmin, vmax = robust_limits(data, symmetric=symmetric)
        for axis in range(3):
            ax = axs[row_idx, axis]
            ax.imshow(slice2d(data, axis, slices[axis]), cmap=cmap, vmin=vmin, vmax=vmax, origin="lower")
            if contour is not None:
                c = slice2d(contour.astype(float), axis, slices[axis])
                if np.any(c > 0):
                    ax.contour(c, levels=[0.5], colors=[contour_color], linewidths=0.8)
            if row_idx == 0:
                ax.set_title(axes_names[axis], fontsize=10)
            if axis == 0:
                ax.set_ylabel(label, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def draw_three_slice_png(
    data: np.ndarray,
    out_png: Path,
    title: str,
    contour: np.ndarray | None = None,
    slices: tuple[int, int, int] | None = None,
    label: str | None = None,
    figsize: tuple[float, float] = (10.5, 3.4),
) -> None:
    if slices is None:
        slices = center_slices(data, contour)
    fig, axs = plt.subplots(1, 3, figsize=figsize, facecolor="white", constrained_layout=True)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    vmin, vmax = robust_limits(data)
    for axis, ax in enumerate(axs):
        ax.imshow(slice2d(data, axis, slices[axis]), cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
        if contour is not None:
            c = slice2d(contour.astype(float), axis, slices[axis])
            if np.nanmin(c) < 0.5 < np.nanmax(c):
                ax.contour(c, levels=[0.5], colors=["#d62728"], linewidths=0.85)
        ax.set_title(["Sagittal", "Coronal", "Axial"][axis], fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    if label:
        fig.text(
            0.015,
            0.94,
            label,
            color="white",
            fontsize=12,
            fontweight="bold",
            bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.72, "pad": 5},
        )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def image_from_three_slice(
    data: np.ndarray,
    title: str,
    contour: np.ndarray | None,
    slices: tuple[int, int, int],
    label: str,
) -> Image.Image:
    import io

    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.4), facecolor="white", constrained_layout=True)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    vmin, vmax = robust_limits(data)
    for axis, ax in enumerate(axs):
        ax.imshow(slice2d(data, axis, slices[axis]), cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
        if contour is not None:
            c = slice2d(contour.astype(float), axis, slices[axis])
            if np.nanmin(c) < 0.5 < np.nanmax(c):
                ax.contour(c, levels=[0.5], colors=["#d62728"], linewidths=0.85)
        ax.set_title(["Sagittal", "Coronal", "Axial"][axis], fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.text(
        0.015,
        0.94,
        label,
        color="white",
        fontsize=12,
        fontweight="bold",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.72, "pad": 5},
    )
    bio = io.BytesIO()
    fig.savefig(bio, format="png", dpi=150)
    plt.close(fig)
    bio.seek(0)
    return Image.open(bio).convert("RGB")


def fade_frames(images: list[Image.Image], hold: int = 2, steps: int = 7) -> tuple[list[Image.Image], list[int]]:
    frames: list[Image.Image] = []
    durations: list[int] = []
    for idx, image in enumerate(images):
        for _ in range(hold):
            frames.append(image)
            durations.append(550 if idx else 850)
        nxt = images[(idx + 1) % len(images)]
        for step in range(1, steps + 1):
            alpha = step / (steps + 1)
            frames.append(Image.blend(image, nxt, alpha))
            durations.append(120)
    return frames, durations


def make_flicker_gif(before_png: Path, after_png: Path, out_gif: Path) -> None:
    before = Image.open(before_png).convert("RGB")
    after = Image.open(after_png).convert("RGB")
    if before.size != after.size:
        after = after.resize(before.size, Image.Resampling.BICUBIC)
    frames = [before, before, after, after]
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_gif,
        save_all=True,
        append_images=frames[1:],
        duration=[650, 250, 650, 250],
        loop=0,
        optimize=True,
    )


def make_flicker_gif_from_images(before: Image.Image, after: Image.Image, out_gif: Path) -> None:
    before = before.convert("RGB")
    after = after.convert("RGB")
    if before.size != after.size:
        after = after.resize(before.size, Image.Resampling.BICUBIC)
    frames = [before, before, after, after]
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_gif,
        save_all=True,
        append_images=frames[1:],
        duration=[700, 250, 700, 250],
        loop=0,
        optimize=True,
    )


def save_single_row_png(
    data: np.ndarray,
    out_png: Path,
    title: str,
    contour: np.ndarray | None = None,
    cmap: str = "gray",
    slices: tuple[int, int, int] | None = None,
) -> None:
    plot_rows([(title, data, cmap, False)], out_png, title, contour=contour, slices=slices)


def corrcoef(a: np.ndarray, b: np.ndarray, mask: np.ndarray | None = None) -> float:
    if mask is None:
        mask = np.isfinite(a) & np.isfinite(b)
    x = a[mask].astype(np.float64)
    y = b[mask].astype(np.float64)
    if x.size < 10:
        return float("nan")
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.sum(x * x) * np.sum(y * y))
    if denom <= 0:
        return float("nan")
    return float(np.sum(x * y) / denom)


def pairwise_corr_matrix(items: list[np.ndarray], mask: np.ndarray | None = None) -> np.ndarray:
    n = len(items)
    sim = np.full((n, n), np.nan, dtype=float)
    for i, data_i in enumerate(items):
        for j, data_j in enumerate(items):
            sim[i, j] = corrcoef(data_i, data_j, mask=mask)
    return sim


def mean_similarity_to_others(sim: np.ndarray) -> list[float]:
    values: list[float] = []
    for idx in range(sim.shape[0]):
        row = np.delete(sim[idx, :], idx)
        values.append(float(np.nanmean(row)) if np.any(np.isfinite(row)) else float("nan"))
    return values


def write_run_sdc(
    run: RunInfo,
    reportlets: Path,
    subject_dir: Path,
    func_name: str,
    xfms_name: str,
    anatomy_ref: Path | None = None,
    contour: np.ndarray | None = None,
    anatomy_slices: tuple[int, int, int] | None = None,
) -> dict[str, str | float | int]:
    if run.distorted is None:
        raise ValueError(f"Run {run.tag} has no distorted/source SBRef")
    distorted_native = load_3d(run.distorted)
    unwarped_native = load_3d(run.unwarped) if run.unwarped is not None else (
        load_3d(run.corrected_acpc) if run.corrected_acpc is not None else distorted_native
    )
    brain = run.run_dir / "MEDIC" / "reference" / "SBref_unwarped_brain.nii.gz"
    if brain.exists():
        mask = load_3d(brain) > 1e-6
    else:
        mask = make_mask(unwarped_native)

    display_contour = None
    if anatomy_ref is not None and contour is not None and anatomy_slices is not None:
        try:
            if run.method == "medic":
                resampled = resample_medic_sdc_pair_to_acpc(run, anatomy_ref, reportlets / "nifti")
            else:
                resampled = resample_standard_sdc_pair_to_acpc(
                    run, subject_dir, func_name, xfms_name, anatomy_ref, reportlets / "nifti"
                )
        except subprocess.CalledProcessError:
            resampled = None
        if resampled is not None:
            distorted = load_3d(resampled[0])
            unwarped = load_3d(resampled[1])
            slices = anatomy_slices
            display_contour = contour
        else:
            distorted = distorted_native
            unwarped = unwarped_native
            slices = center_slices(unwarped_native, mask)
    else:
        distorted = distorted_native
        unwarped = unwarped_native
        slices = center_slices(unwarped_native, mask)

    display_mask = display_contour if display_contour is not None and display_contour.shape == unwarped.shape else make_mask(unwarped)
    distorted_unwarped_r = corrcoef(distorted, unwarped, mask=display_mask) if distorted.shape == unwarped.shape else float("nan")

    base = f"subj_desc-sdc_{run.tag}"
    gif = reportlets / f"{base}_flicker.gif"
    before = image_from_three_slice(
        distorted,
        "Fieldmap correction: distorted vs unwarped",
        display_contour,
        slices,
        f"S{run.session} R{run.run} distorted",
    )
    after = image_from_three_slice(
        unwarped,
        "Fieldmap correction: distorted vs unwarped",
        display_contour,
        slices,
        f"S{run.session} R{run.run} SDC corrected",
    )
    make_flicker_gif_from_images(before, after, gif)

    disp_summary: dict[str, float] = {}
    if (run.run_dir / "MEDIC" / "reference" / "SBref_medic_reference_displacement.nii.gz").exists():
        disp = load_3d(run.run_dir / "MEDIC" / "reference" / "SBref_medic_reference_displacement.nii.gz")
        disp_summary = {
            "disp_p95": float(np.percentile(np.abs(disp[np.isfinite(disp)]), 95)),
            "disp_max": float(np.max(np.abs(disp[np.isfinite(disp)]))),
        }
    elif run.displacement is not None and run.displacement.exists():
        disp = load_3d(run.displacement)
        finite = disp[np.isfinite(disp)]
        if finite.size:
            disp_summary = {
                "disp_p95": float(np.percentile(np.abs(finite), 95)),
                "disp_max": float(np.max(np.abs(finite))),
            }

    return {
        "session": run.session,
        "run": run.run,
        "method": run.method,
        "pe_axis": run.pe_axis,
        "sdc_gif": gif.name,
        "distorted_unwarped_r": distorted_unwarped_r,
        **disp_summary,
    }


def write_longitudinal_mosaic(runs: list[RunInfo], reportlets: Path) -> dict[str, str]:
    images = [(run, load_3d(run.unwarped)) for run in runs]
    masks = [make_mask(data) for _, data in images]
    center = center_slices(images[0][1], masks[0])
    axis = 2
    idx = center[axis]
    ncols = min(4, len(images))
    nrows = int(np.ceil(len(images) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), facecolor="white")
    axs = np.atleast_1d(axs).ravel()
    for ax, (run, data), _mask in zip(axs, images, masks):
        vmin, vmax = robust_limits(data)
        ax.imshow(slice2d(data, axis, idx), cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
        ax.set_title(run.tag, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axs[len(images) :]:
        ax.axis("off")
    fig.suptitle("Longitudinal unwarped SBRef consistency", fontsize=14, fontweight="bold")
    out = reportlets / "subj_desc-longitudinalUnwarped_mosaic.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

    sim = pairwise_corr_matrix([data for _run, data in images])
    mean_to_others = mean_similarity_to_others(sim)
    rows = []
    for (run, _data), mean_r in zip(images, mean_to_others):
        rows.append({"session": run.session, "run": run.run, "mean_unwarped_r_to_other_runs": mean_r})
    return {"mosaic": out.name, "metrics": rows, "run_similarity": sim.tolist()}


def write_run_similarity_matrix(sim: np.ndarray, labels: list[str], out: Path) -> None:
    """Render run similarity as a color-only heatmap, with run labels and scale."""
    n = len(labels)
    fig, ax = plt.subplots(figsize=(5.2, 4.5), facecolor="white", constrained_layout=True)
    im = ax.imshow(sim, cmap="viridis", vmin=np.nanmin(sim), vmax=1.0)
    ax.set_xticks(range(n), labels=labels, rotation=45, ha="right")
    ax.set_yticks(range(n), labels=labels)
    ax.set_title("Run-to-run unwarped SBRef similarity\nafter coregistration")
    # Numeric cell annotations obscure the pattern for multi-session datasets.
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.savefig(out, dpi=150)
    plt.close(fig)


def write_coreg_animation_and_stability(
    subject_dir: Path,
    func_name: str,
    xfms_name: str,
    func_prefix: str,
    runs: list[RunInfo],
    reportlets: Path,
) -> dict[str, object]:
    xfm_dir = subject_dir / "func" / "xfms" / xfms_name
    t1 = xfm_dir / "T1w_acpc_func.nii.gz"
    t2 = xfm_dir / "T2w_acpc_func.nii.gz"
    ribbon = xfm_dir / "CorticalRibbon_acpc_func_mask.nii.gz"
    if not t1.exists() or not ribbon.exists():
        return {}

    t1d = load_3d(t1)
    t2d = load_3d(t2) if t2.exists() else None
    contour = load_3d(ribbon) > 0
    slices = off_midline_sagittal_slices(t1d, contour)
    mask = contour | make_mask(t1d)

    anchor_png = reportlets / "subj_desc-coreg_T1wAnchor.png"
    draw_three_slice_png(t1d, anchor_png, "T1w anatomy anchor", contour=contour, slices=slices, label="T1w anchor")

    sbref_items: list[tuple[RunInfo, np.ndarray, Path]] = []
    resampled_dir = reportlets / "nifti"
    for run in runs:
        sbref_path = resample_run_unwarped_sbref_to_acpc(run, t1, resampled_dir)
        if sbref_path is None:
            continue
        sbref_items.append((run, load_3d(sbref_path), sbref_path))

    frames = [image_from_three_slice(t1d, "Coregistration over time", contour, slices, "T1w anchor")]
    frame_pngs: list[str] = [anchor_png.name]
    for run, sbref, _path in sbref_items:
        label = f"S{run.session} R{run.run} SDC SBRef -> avg + BBR"
        png = reportlets / f"subj_desc-coreg_{run.tag}_unwarpedSBRefBBR.png"
        draw_three_slice_png(sbref, png, "Coregistration over time", contour=contour, slices=slices, label=label)
        frame_pngs.append(png.name)
        frames.append(image_from_three_slice(sbref, "Coregistration over time", contour, slices, label))

    gif = reportlets / "subj_desc-coreg_unwarpedSBRefBBR_fade.gif"
    if frames:
        gif_frames, durations = fade_frames(frames)
        gif_frames[0].save(
            gif,
            save_all=True,
            append_images=gif_frames[1:],
            duration=durations,
            loop=0,
            optimize=True,
        )

    metrics: list[dict[str, object]] = []
    for run, sbref, path in sbref_items:
        row: dict[str, object] = {
            "session": run.session,
            "run": run.run,
            "unwarped_sbref_acpc_file": str(path),
            "corr_sbref_to_t1w": corrcoef(sbref, t1d, mask=mask),
        }
        if t2d is not None:
            row["corr_sbref_to_t2w"] = corrcoef(sbref, t2d, mask=mask)
        metrics.append(row)

    n = len(sbref_items)
    if n == 0:
        return {}
    sim = pairwise_corr_matrix([sbref for _run, sbref, _path in sbref_items], mask=mask)
    mean_sbref_to_others = mean_similarity_to_others(sim)
    for row, mean_r in zip(metrics, mean_sbref_to_others):
        row["mean_sbref_r_to_other_runs"] = mean_r

    labels = [f"S{run.session}R{run.run}" for run, _sbref, _path in sbref_items]
    sim_png = reportlets / "subj_desc-coreg_runSimilarityMatrix.png"
    write_run_similarity_matrix(sim, labels, sim_png)

    stability_png = reportlets / "subj_desc-coreg_anatSimilarityOverTime.png"
    x = np.arange(1, len(metrics) + 1)
    fig, (ax, ax2) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        facecolor="white",
        constrained_layout=True,
        sharex=True,
    )
    if metrics:
        ax.plot(x, [float(m["corr_sbref_to_t1w"]) for m in metrics], marker="o", label="Pearson r to T1w")
        if "corr_sbref_to_t2w" in metrics[0]:
            ax.plot(x, [float(m["corr_sbref_to_t2w"]) for m in metrics], marker="o", label="Pearson r to T2w")
        ax2.plot(
            x,
            [float(m["mean_sbref_r_to_other_runs"]) for m in metrics],
            marker="o",
            color="#2ca02c",
            label="Mean SBRef r to other runs",
        )
    ax.set_xticks(x, labels=[f"S{m['session']}R{m['run']}" for m in metrics])
    ax.set_ylabel("Pearson r")
    ax.set_title("Unwarped SBRef-to-anatomy similarity after coregistration")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    ax2.set_ylabel("Mean pairwise r")
    ax2.set_xlabel("Scan order")
    ax2.set_title("Run reference consistency after coregistration")
    ax2.grid(True, axis="y", alpha=0.25)
    ax2.legend(loc="best")
    fig.savefig(stability_png, dpi=150)
    plt.close(fig)

    return {
        "anchor_png": anchor_png.name,
        "fade_gif": gif.name,
        "frame_pngs": frame_pngs,
        "similarity_matrix_png": sim_png.name,
        "stability_png": stability_png.name,
        "metrics": metrics,
        "run_similarity": sim.tolist(),
    }


def write_coreg(subject_dir: Path, xfms_name: str, reportlets: Path) -> dict[str, str] | None:
    xfm_dir = subject_dir / "func" / "xfms" / xfms_name
    t1 = xfm_dir / "T1w_acpc_func.nii.gz"
    medic_epi = xfm_dir / "AvgMEDICSBref2acpc_EpiReg+BBR.nii.gz"
    standard_epi = xfm_dir / "AvgSBref2acpc_EpiReg+BBR.nii.gz"
    epi = medic_epi if medic_epi.exists() else standard_epi
    ribbon = xfm_dir / "CorticalRibbon_acpc_func_mask.nii.gz"
    if not t1.exists() or not epi.exists():
        return None
    t1d = load_3d(t1)
    epid = load_3d(epi)
    rows = [("T1w functional grid", t1d, "gray", False)]
    rows.append(("EPI after BBR", epid, "gray", False))
    contour = load_3d(ribbon) > 0 if ribbon.exists() else None
    out = reportlets / "subj_desc-coreg_avgReference.png"
    slices = off_midline_sagittal_slices(t1d, contour) if contour is not None else None
    plot_rows(rows, out, "Functional to structural coregistration", contour=contour, slices=slices)
    return {"coreg_panel": out.name}


def write_metrics_csv(metrics: list[dict[str, object]], out_csv: Path) -> None:
    keys: list[str] = []
    for row in metrics:
        for key in row:
            if key not in keys:
                keys.append(key)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    delimiter = "\t" if out_csv.suffix == ".tsv" else ","
    with out_csv.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=keys, extrasaction="ignore", delimiter=delimiter)
        writer.writeheader()
        for row in metrics:
            writer.writerow(row)


def cleanup_obsolete_outputs(coreg_dir: Path, fieldmap_dir: Path) -> None:
    obsolete_files = [
        coreg_dir / "longitudinal_metrics.tsv",
        fieldmap_dir / "longitudinal_metrics.tsv",
        coreg_dir / "reportlets" / "subj_desc-longitudinalUnwarped_mosaic.png",
        fieldmap_dir / "reportlets" / "subj_desc-longitudinalUnwarped_mosaic.png",
        coreg_dir / "reportlets" / "subj_desc-coreg_avgMEDICSBref_legacyRows.png",
        coreg_dir / "index.html",
        fieldmap_dir / "index.html",
    ]
    for path in obsolete_files:
        try:
            path.unlink()
        except FileNotFoundError:
            pass

    fieldmap_reportlets = fieldmap_dir / "reportlets"
    if fieldmap_reportlets.exists():
        for path in fieldmap_reportlets.glob("subj_desc-sdc_*.png"):
            path.unlink()


def cleanup_legacy_experimental_dirs(coreg_dir: Path, fieldmap_dir: Path) -> None:
    for path in (coreg_dir / "FMRIPrepLikeExperimental", fieldmap_dir / "FMRIPrepLikeExperimental"):
        if not path.exists():
            continue
        if path.is_dir():
            shutil.rmtree(path)


def rel_link(from_dir: Path, to_path: Path) -> str:
    return Path(os.path.relpath(to_path, start=from_dir)).as_posix()


def write_coreg_html(
    subject_dir: Path,
    func_name: str,
    coreg_dir: Path,
    fieldmap_report: Path,
    coreg: dict[str, str] | None,
    coreg_time: dict[str, object],
) -> None:
    coreg_dir.mkdir(parents=True, exist_ok=True)
    coreg_html = ""
    if coreg:
        coreg_html = (
            "<section><h2>Functional-Anatomical Coregistration: Average Reference</h2>"
            f"<img src='reportlets/{html.escape(coreg['coreg_panel'])}' alt='Coregistration panel'>"
            "</section>"
        )

    gif_html = ""
    matrix_html = ""
    pearson_html = ""
    if coreg_time:
        gif_html = (
            "<section><h2>Runwise Unwarped SBRef+BBR Alignment to Anatomy</h2>"
            f"<img src='reportlets/{html.escape(str(coreg_time['fade_gif']))}' alt='Runwise unwarped SBRef BBR fade animation'>"
            "</section>"
        )
        matrix_html = (
            "<section><h2>Run-to-Run Similarity Matrix</h2>"
            f"<img src='reportlets/{html.escape(str(coreg_time['similarity_matrix_png']))}' alt='Run-to-run unwarped SBRef similarity matrix'>"
            "</section>"
        )
        pearson_html = (
            "<section><h2>Pearson Registration Stability</h2>"
            f"<img src='reportlets/{html.escape(str(coreg_time['stability_png']))}' alt='Unwarped SBRef Pearson similarity over time'>"
            "</section>"
        )

    fieldmap_href = rel_link(coreg_dir, fieldmap_report)
    report = coreg_dir / COREG_REPORT_FILENAME
    report.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(subject_dir.name)} {html.escape(func_name)} Coreg QC</title>
<style>
body {{ margin: 0; font-family: Arial, sans-serif; color: #1f2328; background: #f6f7f9; }}
header {{ padding: 22px 28px; background: #ffffff; border-bottom: 1px solid #d8dee4; }}
main {{ max-width: 1180px; margin: 0 auto; padding: 22px; }}
section {{ background: #ffffff; border: 1px solid #d8dee4; border-radius: 6px; padding: 16px; margin-bottom: 18px; }}
h1, h2 {{ margin: 0 0 12px; }}
img {{ max-width: 100%; height: auto; display: block; margin: 10px 0; border: 1px solid #d8dee4; }}
a.button {{ display: inline-block; padding: 8px 11px; border: 1px solid #8c959f; border-radius: 5px; color: #1f2328; text-decoration: none; background: #fff; }}
section img + img {{ margin-top: 16px; }}
</style>
</head>
<body>
<header>
<h1>{html.escape(subject_dir.name)} {html.escape(func_name)} Coreg QC</h1>
<div>Experimental fMRIPrep-style coregistration visuals.</div>
</header>
<main>
<section><a class="button" href="{html.escape(fieldmap_href)}">Open Fieldmap/SDC Report</a></section>
{coreg_html}
{gif_html}
{matrix_html}
{pearson_html}
</main>
</body>
</html>
""",
        encoding="utf-8",
    )


def write_fieldmap_html(
    subject_dir: Path,
    func_name: str,
    fieldmap_dir: Path,
    coreg_report: Path,
    run_metrics: list[dict[str, object]],
) -> None:
    rows = []
    for row in run_metrics:
        rows.append(
            "<tr>"
            f"<td>S{row['session']} R{row['run']}</td>"
            f"<td>{html.escape(str(row.get('method', '')))}</td>"
            f"<td>{html.escape(str(row.get('pe_axis', '')))}</td>"
            f"<td>{float(row.get('distorted_unwarped_r', float('nan'))):.3f}</td>"
            f"<td>{float(row.get('disp_p95', float('nan'))):.3f}</td>"
            f"<td><a href='reportlets/{html.escape(str(row['sdc_gif']))}'>flicker</a></td>"
            "</tr>"
        )

    run_cards = []
    for row in run_metrics:
        run_cards.append(
            "<section class='run'>"
            f"<h3>S{row['session']} R{row['run']}</h3>"
            f"<img src='reportlets/{html.escape(str(row['sdc_gif']))}' alt='SDC flicker {row['session']} {row['run']}'>"
            "</section>"
        )

    report = fieldmap_dir / FIELDMAP_REPORT_FILENAME
    coreg_href = rel_link(fieldmap_dir, coreg_report)
    report.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(subject_dir.name)} {html.escape(func_name)} Fieldmap/SDC QC</title>
<style>
body {{ margin: 0; font-family: Arial, sans-serif; color: #1f2328; background: #f6f7f9; }}
header {{ padding: 22px 28px; background: #ffffff; border-bottom: 1px solid #d8dee4; }}
main {{ max-width: 1180px; margin: 0 auto; padding: 22px; }}
section {{ background: #ffffff; border: 1px solid #d8dee4; border-radius: 6px; padding: 16px; margin-bottom: 18px; }}
h1, h2, h3 {{ margin: 0 0 12px; }}
img {{ max-width: 100%; height: auto; display: block; margin: 10px 0; border: 1px solid #d8dee4; }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ border-bottom: 1px solid #d8dee4; padding: 8px 10px; text-align: left; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 18px; }}
.run img:first-of-type {{ width: 100%; }}
code {{ background: #eef1f4; padding: 1px 4px; border-radius: 4px; }}
a.button {{ display: inline-block; padding: 8px 11px; border: 1px solid #8c959f; border-radius: 5px; color: #1f2328; text-decoration: none; background: #fff; }}
</style>
</head>
<body>
<header>
<h1>{html.escape(subject_dir.name)} {html.escape(func_name)} Fieldmap/SDC QC</h1>
<div>Experimental fMRIPrep-style distortion-correction visuals from existing pipeline outputs.</div>
</header>
<main>
<section><a class="button" href="{html.escape(coreg_href)}">Open Coreg Report</a></section>
<section>
<h2>Run Summary</h2>
<table>
<thead><tr><th>Run</th><th>Method</th><th>PE axis</th><th>Before/After r</th><th>|field/disp| p95</th><th>Flicker</th></tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>
</section>
<div class="grid">
{''.join(run_cards)}
</div>
</main>
</body>
</html>
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject-dir", required=True, type=Path)
    parser.add_argument("--func-name", default="rest")
    parser.add_argument("--xfms-name", help="Transform namespace under func/xfms; defaults to --func-name")
    parser.add_argument("--func-prefix", default="Rest")
    parser.add_argument("--distortion-mode", default="auto", choices=["auto", "topup", "direct_b0", "fieldmap", "gre_fieldmap", "medic", "none"])
    parser.add_argument("--out-dir", type=Path, help="Write one combined report tree instead of split Coreg/Fieldmap trees")
    parser.add_argument("--coreg-dir", type=Path)
    parser.add_argument("--fieldmap-dir", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    subject_dir = args.subject_dir.resolve()
    xfms_name = args.xfms_name or args.func_name
    if args.out_dir:
        coreg_dir = args.out_dir.resolve()
        fieldmap_dir = args.out_dir.resolve()
    else:
        qc_dirname = os.environ.get("QC_OUTPUT_DIRNAME", "QC")
        qa_dir = subject_dir / "func" / args.func_name / qc_dirname
        coreg_dir = (args.coreg_dir or qa_dir / "Coreg").resolve()
        fieldmap_dir = (args.fieldmap_dir or qa_dir / "Fieldmaps").resolve()
    coreg_reportlets = coreg_dir / "reportlets"
    fieldmap_reportlets = fieldmap_dir / "reportlets"
    coreg_reportlets.mkdir(parents=True, exist_ok=True)
    fieldmap_reportlets.mkdir(parents=True, exist_ok=True)
    cleanup_obsolete_outputs(coreg_dir, fieldmap_dir)
    cleanup_legacy_experimental_dirs(coreg_dir, fieldmap_dir)

    runs = discover_runs(subject_dir, args.func_name, xfms_name, args.distortion_mode)
    if not runs:
        raise SystemExit(f"No SDC/coregistration run-level SBRef outputs found under {subject_dir}/func/{args.func_name}")

    xfm_dir = subject_dir / "func" / "xfms" / xfms_name
    fieldmap_anatomy_ref = xfm_dir / "T1w_acpc_func.nii.gz"
    fieldmap_ribbon = xfm_dir / "CorticalRibbon_acpc_func_mask.nii.gz"
    fieldmap_contour = None
    fieldmap_slices = None
    if fieldmap_anatomy_ref.exists() and fieldmap_ribbon.exists():
        anatomy = load_3d(fieldmap_anatomy_ref)
        fieldmap_contour = load_3d(fieldmap_ribbon) > 0
        fieldmap_slices = center_slices(anatomy, fieldmap_contour)
    else:
        fieldmap_anatomy_ref = None

    run_metrics = [
        write_run_sdc(
            run,
            fieldmap_reportlets,
            subject_dir,
            args.func_name,
            xfms_name,
            fieldmap_anatomy_ref,
            fieldmap_contour,
            fieldmap_slices,
        )
        for run in runs
    ]
    coreg = write_coreg(subject_dir, xfms_name, coreg_reportlets)
    coreg_time = write_coreg_animation_and_stability(
        subject_dir,
        args.func_name,
        xfms_name,
        args.func_prefix,
        runs,
        coreg_reportlets,
    )

    payload = {
        "subject_dir": str(subject_dir),
        "func_name": args.func_name,
        "xfms_name": xfms_name,
        "distortion_mode": args.distortion_mode,
        "run_metrics": run_metrics,
        "coreg": coreg,
        "coreg_time": coreg_time,
        "coreg_dir": str(coreg_dir),
        "fieldmap_dir": str(fieldmap_dir),
    }
    (coreg_dir / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    (fieldmap_dir / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    write_metrics_csv(run_metrics, fieldmap_dir / "run_metrics.tsv")
    if coreg_time.get("metrics"):
        write_metrics_csv(coreg_time["metrics"], coreg_dir / "coreg_stability_metrics.tsv")
    write_fieldmap_html(
        subject_dir,
        args.func_name,
        fieldmap_dir,
        coreg_dir / COREG_REPORT_FILENAME,
        run_metrics,
    )
    write_coreg_html(
        subject_dir,
        args.func_name,
        coreg_dir,
        fieldmap_dir / FIELDMAP_REPORT_FILENAME,
        coreg,
        coreg_time,
    )
    print(coreg_dir / COREG_REPORT_FILENAME)
    print(fieldmap_dir / FIELDMAP_REPORT_FILENAME)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

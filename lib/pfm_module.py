#!/usr/bin/env python3
"""PFM preparation module: concatenate run CIFTIs and apply FD censoring.

This mirrors the front end of ExamplePFMCall.m while deferring RidgeFusion/NSI
porting to later stages.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import nibabel as nib
from nibabel.cifti2.cifti2_axes import SeriesAxis


def compute_fd_from_mcfpar(mcf_par: Path, tr: float) -> np.ndarray:
    rp = np.loadtxt(mcf_par, dtype=np.float64)
    if rp.ndim == 1:
        rp = rp.reshape(1, -1)
    if rp.shape[1] != 6:
        raise ValueError(f"Expected 6-column motion file: {mcf_par}, got {rp.shape}")

    # Match legacy MATLAB behavior: backward difference over ~2.5 s.
    n_trs = int(round(2.5 / tr))
    n_trs = max(1, min(n_trs, rp.shape[0]))

    fd = np.zeros_like(rp)
    for col in range(6):
        for i in range(n_trs, rp.shape[0]):
            fd[i, col] = abs(rp[i, col] - rp[i - n_trs, col])

    fd_ang = fd[:, :3] / (2.0 * np.pi)
    fd_ang = fd_ang * 100.0 * np.pi
    fd_lin = fd[:, 3:]
    return np.sum(np.hstack([fd_lin, fd_ang]), axis=1)


def _parse_trailing_int(path_name: str, prefix: str) -> Optional[int]:
    if not path_name.startswith(prefix):
        return None
    try:
        return int(path_name.split("_")[-1])
    except Exception:
        return None


def discover_runs(subdir: Path, func_dirname: str, start_session: int) -> List[Tuple[int, int, Path]]:
    out: List[Tuple[int, int, Path]] = []
    base = subdir / "func" / func_dirname
    sessions = []
    for ses in base.glob("session_*"):
        s = _parse_trailing_int(ses.name, "session_")
        if s is None:
            continue
        sessions.append((s, ses))

    for s, ses in sorted(sessions, key=lambda x: x[0]):
        if s < start_session:
            continue
        runs = []
        for run in ses.glob("run_*"):
            r = _parse_trailing_int(run.name, "run_")
            if r is None:
                continue
            runs.append((r, run))
        for r, run in sorted(runs, key=lambda x: x[0]):
            out.append((s, r, run))
    return out


def resolve_run_cifti(run_dir: Path, func_file_prefix: str, input_tag: str) -> Optional[Path]:
    exact = run_dir / f"{func_file_prefix}_{input_tag}.dtseries.nii"
    if exact.exists():
        return exact

    stamped = sorted(run_dir.glob(f"{func_file_prefix}_{input_tag}_*.dtseries.nii"))
    if not stamped:
        return None
    if len(stamped) == 1:
        return stamped[0]
    raise RuntimeError(
        f"Ambiguous PFM input for {run_dir}: found multiple candidates matching "
        f"{func_file_prefix}_{input_tag}_*.dtseries.nii. "
        "Set --input-tag to an exact tag or remove stale alternatives."
    )


def read_tr(run_dir: Path) -> float:
    tr_txt = run_dir / "TR.txt"
    if not tr_txt.exists():
        return 2.0
    return float(tr_txt.read_text().strip())


def save_dtseries(data_time_by_gray: np.ndarray, ref_img: nib.Cifti2Image, out_path: Path) -> None:
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    if len(axes) != 2:
        raise ValueError("Expected dense time series CIFTI with 2 axes")
    series_axis = axes[0]
    brain_axis = axes[1]
    if not isinstance(series_axis, SeriesAxis):
        raise ValueError("Expected SeriesAxis for dtseries first axis")

    new_series = SeriesAxis(series_axis.start, series_axis.step, data_time_by_gray.shape[0], unit=series_axis.unit)
    new_header = nib.Cifti2Header.from_axes((new_series, brain_axis))
    out_img = nib.Cifti2Image(data_time_by_gray.astype(np.float32), new_header, nifti_header=ref_img.nifti_header)
    nib.save(out_img, str(out_path))


def main() -> int:
    ap = argparse.ArgumentParser(description="PFM prep-only module (concatenate + FD censor)")
    ap.add_argument("--subject", required=True)
    ap.add_argument("--study-folder", required=True)
    ap.add_argument("--mode", default="prep_only")
    ap.add_argument("--start-session", type=int, default=1)
    ap.add_argument("--func-dirname", default="rest")
    ap.add_argument("--func-file-prefix", default="Rest")
    ap.add_argument("--input-tag", default="OCME+MEICA+MGTR")
    ap.add_argument("--fd-threshold", type=float, default=0.3)
    ap.add_argument("--out-subdir", default="ConcatenatedCiftis")
    # PFM feature knobs (many are staged for future full-mode implementation).
    ap.add_argument("--demean-runs", type=int, default=1)
    ap.add_argument("--var-norm-runs", type=int, default=0)
    ap.add_argument("--var-norm-eps", type=float, default=1e-8)
    ap.add_argument("--censor-by-fd", type=int, default=1)
    ap.add_argument("--concatenate-runs", type=int, default=1)
    ap.add_argument("--save-fd-txt", type=int, default=1)
    ap.add_argument("--save-scanidx-txt", type=int, default=1)
    ap.add_argument("--nsi-enable", type=int, default=0)
    ap.add_argument("--nsi-ridge-lambdas", default="1,5,10,25,50")
    ap.add_argument("--nsi-headline-lambda", type=float, default=10.0)
    ap.add_argument("--nsi-compute-morans", type=int, default=0)
    ap.add_argument("--nsi-compute-slope", type=int, default=0)
    ap.add_argument("--rf-enable", type=int, default=0)
    ap.add_argument("--rf-outfile", default="RidgeFusion_VTX")
    ap.add_argument("--rf-fc-weight", type=float, default=1.0)
    ap.add_argument("--rf-spatial-weight", type=float, default=0.1)
    ap.add_argument("--ap-enable", type=int, default=0)
    ap.add_argument("--ap-method", default="kmeans")
    ap.add_argument("--ap-kmax-cap", type=int, default=8)
    ap.add_argument("--ap-kmin-per-patch", type=int, default=1)
    ap.add_argument("--ap-diffuse-iters", type=int, default=6)
    ap.add_argument("--ap-lambda-spatial", type=float, default=8.0)
    ap.add_argument("--priors-mat", default="")
    ap.add_argument("--subcort-priors-nii", default="")
    ap.add_argument("--neighbors-mat", default="")
    args = ap.parse_args()

    subdir = Path(args.study_folder) / args.subject
    runs = discover_runs(subdir, args.func_dirname, args.start_session)
    if not runs:
        raise SystemExit(f"No runs found in {subdir / 'func' / args.func_dirname}")

    outdir = subdir / "func" / args.func_dirname / args.out_subdir
    outdir.mkdir(parents=True, exist_ok=True)

    concat_blocks = []
    fd_all = []
    scan_idx = []
    ref_img = None

    run_meta = []
    missing_runs = []

    for scan_id, (s, r, run_dir) in enumerate(runs, start=1):
        in_dt = resolve_run_cifti(run_dir, args.func_file_prefix, args.input_tag)
        if in_dt is None:
            missing_runs.append(str(run_dir / f"{args.func_file_prefix}_{args.input_tag}.dtseries.nii"))
            continue

        mcf_par = run_dir / "MCF.par"
        if not mcf_par.exists():
            mcf_par = run_dir / "mcf.par"
        if not mcf_par.exists():
            raise FileNotFoundError(f"Missing MCF.par/mcf.par needed for FD: {run_dir}")

        tr = read_tr(run_dir)
        fd = compute_fd_from_mcfpar(mcf_par, tr)

        img = nib.load(str(in_dt))
        data = np.asanyarray(img.dataobj)
        if data.ndim != 2:
            raise ValueError(f"Expected 2D dtseries data (time x grayordinates): {in_dt}, got {data.shape}")
        if data.shape[0] != fd.shape[0]:
            raise ValueError(f"FD length mismatch for {in_dt}: fd={fd.shape[0]} vs time={data.shape[0]}")

        data = data.astype(np.float32, copy=False)
        if int(args.demean_runs) == 1:
            # Mirror MATLAB concatenate_scans behavior: demean each run before concat.
            run_mean = np.mean(data, axis=0, keepdims=True, dtype=np.float64).astype(np.float32)
            data = data - run_mean
        if int(args.var_norm_runs) == 1:
            # Optional run-wise variance normalization after demeaning.
            run_std = np.std(data, axis=0, keepdims=True, dtype=np.float64).astype(np.float32)
            safe_std = np.where(run_std > float(args.var_norm_eps), run_std, 1.0).astype(np.float32)
            data = data / safe_std

        if ref_img is None:
            ref_img = img
        elif data.shape[1] != concat_blocks[0].shape[1]:
            raise ValueError(f"Grayordinate mismatch in {in_dt}")

        concat_blocks.append(data)
        fd_all.append(fd.astype(np.float32, copy=False))
        scan_idx.append(np.full(fd.shape[0], scan_id, dtype=np.int32))
        run_meta.append(
            {
                "scan_id": int(scan_id),
                "session": int(s),
                "run": int(r),
                "run_dir": str(run_dir),
                "input_dtseries": str(in_dt),
                "n_timepoints": int(fd.shape[0]),
            }
        )

    if not concat_blocks:
        raise SystemExit("No valid run dtseries found for PFM prep")
    if missing_runs:
        missing_preview = "\n".join(missing_runs[:12])
        raise SystemExit(
            "PFM prep found run directories with missing CIFTI inputs.\n"
            f"Expected tag: {args.func_file_prefix}_{args.input_tag}.dtseries.nii "
            "(or a single stamped variant with '_<stamp>').\n"
            f"Missing ({len(missing_runs)} runs):\n{missing_preview}"
        )

    if int(args.concatenate_runs) != 1:
        raise SystemExit("PFM prep currently requires --concatenate-runs=1")
    concat = np.vstack(concat_blocks)
    fd_vec = np.concatenate(fd_all)
    scan_vec = np.concatenate(scan_idx)

    if int(args.censor_by_fd) == 1:
        keep = fd_vec < float(args.fd_threshold)
    else:
        keep = np.ones(fd_vec.shape[0], dtype=bool)
    cens = concat[keep, :]

    base = f"{args.func_file_prefix}_{args.input_tag}"
    out_concat = outdir / f"{base}_Concatenated.dtseries.nii"
    out_cens = outdir / f"{base}_Concatenated+FDlt{str(args.fd_threshold).replace('.', 'p')}.dtseries.nii"

    assert ref_img is not None
    save_dtseries(concat, ref_img, out_concat)
    save_dtseries(cens, ref_img, out_cens)

    if int(args.save_fd_txt) == 1:
        np.savetxt(outdir / "FD.txt", fd_vec, fmt="%.6f")
    if int(args.save_scanidx_txt) == 1:
        np.savetxt(outdir / "ScanIdx.txt", scan_vec, fmt="%d")
    (outdir / "ScanInfo.json").write_text(json.dumps(run_meta, indent=2))

    summary = {
        "subject": args.subject,
        "start_session": args.start_session,
        "func_dirname": args.func_dirname,
        "func_file_prefix": args.func_file_prefix,
        "input_tag": args.input_tag,
        "fd_threshold": args.fd_threshold,
        "n_runs_used": int(len(concat_blocks)),
        "n_scans": int(len(run_meta)),
        "n_runs_missing_input": int(len(missing_runs)),
        "n_timepoints_concat": int(concat.shape[0]),
        "n_timepoints_kept": int(cens.shape[0]),
        "frac_kept": float(cens.shape[0] / concat.shape[0]),
        "out_concat": str(out_concat),
        "out_censored": str(out_cens),
        "run_scan_map": run_meta,
        "mode": args.mode,
        "pfm_knobs": {
            "demean_runs": int(args.demean_runs),
            "var_norm_runs": int(args.var_norm_runs),
            "var_norm_eps": float(args.var_norm_eps),
            "censor_by_fd": int(args.censor_by_fd),
            "concatenate_runs": int(args.concatenate_runs),
            "save_fd_txt": int(args.save_fd_txt),
            "save_scanidx_txt": int(args.save_scanidx_txt),
            "nsi_enable": int(args.nsi_enable),
            "nsi_ridge_lambdas": args.nsi_ridge_lambdas,
            "nsi_headline_lambda": float(args.nsi_headline_lambda),
            "nsi_compute_morans": int(args.nsi_compute_morans),
            "nsi_compute_slope": int(args.nsi_compute_slope),
            "rf_enable": int(args.rf_enable),
            "rf_outfile": args.rf_outfile,
            "rf_fc_weight": float(args.rf_fc_weight),
            "rf_spatial_weight": float(args.rf_spatial_weight),
            "ap_enable": int(args.ap_enable),
            "ap_method": args.ap_method,
            "ap_kmax_cap": int(args.ap_kmax_cap),
            "ap_kmin_per_patch": int(args.ap_kmin_per_patch),
            "ap_diffuse_iters": int(args.ap_diffuse_iters),
            "ap_lambda_spatial": float(args.ap_lambda_spatial),
            "priors_mat": args.priors_mat,
            "subcort_priors_nii": args.subcort_priors_nii,
            "neighbors_mat": args.neighbors_mat,
        },
    }
    (outdir / "PFMPrepSummary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

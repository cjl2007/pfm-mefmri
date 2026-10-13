#!/usr/bin/env python3
"""
No-MATLAB replacement for find_epi_params.m.
Creates run-level metadata sidecars used by coreg/headmotion:
  - TE.txt, TR.txt, EffectiveEchoSpacing.txt, SliceTiming.txt, PE.txt
Also updates:
  - func/<FuncName>/QC/Validation/AllTEs.txt, func/<FuncName>/QC/Validation/AllTRs.txt
  - func/xfms/<FuncName>/EffectiveEchoSpacing.txt (when StartSession == 1)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import List


def run(cmd: List[str]) -> str:
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout.strip()


def fmt_num(x: float) -> str:
    # Approximate MATLAB num2str default text output for this pipeline.
    return f"{x:.5g}"


def write_validation_visuals(validation_dir: Path, rows: List[dict]) -> None:
    if not rows:
        return
    try:
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/mefmri_matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except Exception as exc:
        (validation_dir / "ValidationQC_warning.txt").write_text(f"Could not write validation figures: {exc}\n")
        return

    labels = [f"S{r['session']}R{r['run']}E{r['echo']}" for r in rows]
    nvol = [r["nvol"] for r in rows]
    size_mb = [r["size_mb"] for r in rows]
    te_ms = [r["te_ms"] for r in rows]
    tr_s = [r["tr_s"] for r in rows]
    echo_counts = {}
    run_vols = {}
    for r in rows:
        key = (r["session"], r["run"])
        echo_counts[key] = max(echo_counts.get(key, 0), r["echo"])
        run_vols.setdefault(key, []).append(r["nvol"])

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), dpi=150)
    fig.suptitle("Input Validation QC", fontsize=14)
    x = list(range(len(rows)))
    for ax, values, title, ylabel, color in [
        (axes[0, 0], nvol, "Volumes by run/echo", "Volumes", "#4c78a8"),
        (axes[0, 1], size_mb, "File size by run/echo", "MB", "#72b7b2"),
    ]:
        ax.bar(x, values, color=color)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    run_keys = sorted(echo_counts)
    run_labels = [f"S{s}R{r}" for s, r in run_keys]
    run_echo_counts = [echo_counts[k] for k in run_keys]
    axes[1, 0].bar(run_labels, run_echo_counts, color="#f28e2b")
    axes[1, 0].set_title("Echo count by run")
    axes[1, 0].set_ylabel("Echoes")
    axes[1, 0].tick_params(axis="x", rotation=45, labelsize=7)
    axes[1, 0].grid(axis="y", color="#dddddd", linewidth=0.5)

    axes[1, 1].scatter(te_ms, tr_s, s=40, color="#54a24b")
    axes[1, 1].set_title("TE/TR consistency")
    axes[1, 1].set_xlabel("TE (ms)")
    axes[1, 1].set_ylabel("TR (s)")
    axes[1, 1].grid(color="#dddddd", linewidth=0.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(validation_dir / "ValidationSummary.png")
    plt.close(fig)

    flagged = []
    for key, vals in sorted(run_vols.items()):
        if len(set(vals)) > 1:
            flagged.append(f"S{key[0]}R{key[1]} volume mismatch across echoes: {sorted(vals)}")
    max_echo = max(run_echo_counts) if run_echo_counts else 0
    for key, count in sorted(echo_counts.items()):
        if count != max_echo:
            flagged.append(f"S{key[0]}R{key[1]} echo count differs from max: {count} vs {max_echo}")
    status = "No run/echo truncation flags detected." if not flagged else "<br>".join(flagged)
    (validation_dir / "index.html").write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Input Validation QC</title>
<style>body{{font-family:Arial,sans-serif;margin:24px;color:#222}}img{{max-width:100%;border:1px solid #ccc}}.flag{{padding:12px;background:#fff4e6;border:1px solid #f0c27b}}</style>
</head><body><h1>Input Validation QC</h1><div class="flag">{status}</div>
<p>Checks echo count, number of volumes, file size, TE, and TR across all runs/echoes.</p>
<img src="ValidationSummary.png" alt="Input validation summary">
<p>Raw tables: <a href="InputValidation.tsv">InputValidation.tsv</a>, <a href="FileSize.tsv">FileSize.tsv</a>, <a href="NumberOfVolumes.tsv">NumberOfVolumes.tsv</a>, <a href="AllTEs.txt">AllTEs.txt</a>, <a href="AllTRs.txt">AllTRs.txt</a>.</p>
</body></html>
""",
        encoding="utf-8",
    )


def phase_to_pedir(phase: str) -> str:
    # Keep legacy behavior for compatibility.
    return "-y" if phase == "j-" else "y"


def build_slice_timing(slice_times: List[float], mb_factor: int) -> List[float]:
    # Replicate MATLAB logic:
    # 1) Use first N unique slice times from the provided list
    # 2) sort to get slice order
    # 3) choose middle slice as reference (MATLAB round)
    n_unique = len(sorted(set(slice_times)))
    base = slice_times[:n_unique]
    # 1-based order indices like MATLAB sort output.
    order_pairs = sorted(enumerate(base, start=1), key=lambda p: p[1])
    slice_order = [idx for idx, _ in order_pairs]
    n = len(slice_order)
    ref_pos = int(math.floor(n / 2 + 0.5))
    ref_slice = slice_order[ref_pos - 1]

    timing = []
    for i in range(1, n + 1):
        pos_ref = slice_order.index(ref_slice) + 1
        pos_i = slice_order.index(i) + 1
        timing.append((pos_ref - pos_i) / n)

    if mb_factor > 1:
        timing = timing * mb_factor
    return timing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subdir", required=True)
    ap.add_argument("--func-name", default="rest")
    ap.add_argument("--func-prefix", default="Rest")
    ap.add_argument("--start-session", type=int, default=1)
    ap.add_argument("--no-fieldmap-mode", action="store_true")
    ap.add_argument("--qc-dirname", default=os.environ.get("QC_OUTPUT_DIRNAME", "QC"))
    ap.add_argument("--validation-qc", type=int, default=int(os.environ.get("QC_VALIDATION_ENABLE", "1")))
    args = ap.parse_args()

    subdir = Path(args.subdir)
    func_name = args.func_name
    func_prefix = args.func_prefix
    start_session = args.start_session
    in_root = subdir / "func" / "unprocessed" / func_name
    out_root = subdir / "func" / func_name
    qa_dir = subdir / "func" / func_name / args.qc_dirname / "Validation"
    xfms_dir = subdir / "func" / "xfms" / func_name

    if start_session == 1 and out_root.exists():
        # Do not delete the whole task directory: fieldmaps now live under
        # func/<task>/field_maps and are produced before coreg calls this helper.
        for stale in out_root.glob("session_*"):
            if stale.is_dir():
                shutil.rmtree(stale)
        avg_sbref = out_root / "AverageSBref"
        if avg_sbref.exists():
            shutil.rmtree(avg_sbref)

    sessions = sorted([p for p in in_root.glob("session_*") if p.is_dir()], key=lambda p: int(p.name.split("_")[1]))
    all_ees = []
    te_lines = []
    tr_lines = []

    qa_dir.mkdir(parents=True, exist_ok=True)
    xfms_dir.mkdir(parents=True, exist_ok=True)

    for sdir in sessions:
        s_idx = int(sdir.name.split("_")[1])
        if s_idx < start_session:
            continue
        runs = sorted([p for p in sdir.glob("run_*") if p.is_dir()], key=lambda p: int(p.name.split("_")[1]))
        for rdir in runs:
            r_idx = int(rdir.name.split("_")[1])
            out_run = out_root / f"session_{s_idx}" / f"run_{r_idx}"
            out_run.mkdir(parents=True, exist_ok=True)

            mag_json_re = re.compile(rf"^{re.escape(func_prefix)}_S\d+_R\d+_E\d+\.json$")
            json_files = sorted([p for p in rdir.glob(f"{func_prefix}*.json") if mag_json_re.match(p.name)])
            if not json_files:
                raise FileNotFoundError(f"No {func_prefix}*.json files in {rdir}")

            te_vals = []
            last_j = None
            for jf in json_files:
                j = json.loads(jf.read_text())
                last_j = j
                te_vals.append(float(j["EchoTime"]) * 1e3)

            assert last_j is not None
            tr = float(last_j["RepetitionTime"])
            ees_val = last_j.get("EffectiveEchoSpacing")
            ees = None if ees_val is None else float(ees_val)
            if ees is not None:
                all_ees.append(ees)

            first_j = json.loads(json_files[0].read_text())
            slice_times = [float(x) for x in first_j["SliceTiming"]]
            mb = int(first_j.get("MultibandAccelerationFactor", 1))
            slice_timing = build_slice_timing(slice_times, mb)
            phase = str(last_j.get("PhaseEncodingDirection", ""))
            pedir = phase_to_pedir(phase)

            (out_run / "TE.txt").write_text(" ".join(fmt_num(v) for v in te_vals) + "\n")
            (out_run / "TR.txt").write_text(fmt_num(tr) + "\n")
            if ees is not None:
                (out_run / "EffectiveEchoSpacing.txt").write_text(fmt_num(ees) + "\n")
            elif not args.no_fieldmap_mode:
                raise KeyError("EffectiveEchoSpacing")
            (out_run / "SliceTiming.txt").write_text("\n".join(fmt_num(v) for v in slice_timing) + "\n")
            (out_run / "PE.txt").write_text(pedir + "\n")

            te_lines.append(f"Session {s_idx} Run {r_idx}: " + " ".join(fmt_num(v) for v in te_vals))
            tr_lines.append(f"Session {s_idx} Run {r_idx}: {fmt_num(tr)}")

    mode_ees = Counter(all_ees).most_common(1)
    if mode_ees and start_session == 1:
        (xfms_dir / "EffectiveEchoSpacing.txt").write_text(fmt_num(mode_ees[0][0]) + "\n")

    all_te_file = qa_dir / "AllTEs.txt"
    all_tr_file = qa_dir / "AllTRs.txt"
    te_text = "\n".join(te_lines) + ("\n" if te_lines else "")
    tr_text = "\n".join(tr_lines) + ("\n" if tr_lines else "")
    if start_session == 1:
        all_te_file.write_text(te_text)
        all_tr_file.write_text(tr_text)
    else:
        with all_te_file.open("a") as f:
            f.write(te_text)
        with all_tr_file.open("a") as f:
            f.write(tr_text)

    # Optional QA matrices replacement: write simple TSVs (MATLAB .mat not required by coreg/headmotion).
    # Keep light-weight and dependency-free.
    fsize_tsv = qa_dir / "FileSize.tsv"
    nvol_tsv = qa_dir / "NumberOfVolumes.tsv"
    rows_fsize = []
    rows_nvol = []
    rows_validation = []
    for sdir in sessions:
        s_idx = int(sdir.name.split("_")[1])
        if s_idx < start_session:
            continue
        runs = sorted([p for p in sdir.glob("run_*") if p.is_dir()], key=lambda p: int(p.name.split("_")[1]))
        for rdir in runs:
            r_idx = int(rdir.name.split("_")[1])
            mag_nii_re = re.compile(rf"^{re.escape(func_prefix)}_S\d+_R\d+_E\d+\.nii\.gz$")
            nii_files = sorted([p for p in rdir.glob(f"{func_prefix}_S*_R*_E*.nii.gz") if mag_nii_re.match(p.name)])
            for e_idx, nf in enumerate(nii_files, start=1):
                size_mb = nf.stat().st_size / 1e6
                nvol = int(run(["fslnvols", str(nf)]))
                rows_fsize.append(f"{s_idx}\t{r_idx}\t{e_idx}\t{size_mb:.6f}")
                rows_nvol.append(f"{s_idx}\t{r_idx}\t{e_idx}\t{nvol}")
                te_ms = math.nan
                tr_s = math.nan
                jf = nf.with_suffix("").with_suffix(".json")
                if jf.exists():
                    j = json.loads(jf.read_text())
                    te_ms = float(j.get("EchoTime", math.nan)) * 1e3
                    tr_s = float(j.get("RepetitionTime", math.nan))
                rows_validation.append(
                    {
                        "session": s_idx,
                        "run": r_idx,
                        "echo": e_idx,
                        "size_mb": size_mb,
                        "nvol": nvol,
                        "te_ms": te_ms,
                        "tr_s": tr_s,
                        "file": nf.name,
                    }
                )
    fsize_tsv.write_text("session\trun\techo\tsize_mb\n" + "\n".join(rows_fsize) + ("\n" if rows_fsize else ""))
    nvol_tsv.write_text("session\trun\techo\tnvol\n" + "\n".join(rows_nvol) + ("\n" if rows_nvol else ""))
    validation_lines = [
        f"{r['session']}\t{r['run']}\t{r['echo']}\t{r['nvol']}\t{r['size_mb']:.6f}\t{fmt_num(r['te_ms'])}\t{fmt_num(r['tr_s'])}\t{r['file']}"
        for r in rows_validation
    ]
    (qa_dir / "InputValidation.tsv").write_text(
        "session\trun\techo\tnvol\tsize_mb\tte_ms\ttr_s\tfile\n"
        + "\n".join(validation_lines)
        + ("\n" if validation_lines else "")
    )
    if int(args.validation_qc) == 1:
        write_validation_visuals(qa_dir, rows_validation)


if __name__ == "__main__":
    main()

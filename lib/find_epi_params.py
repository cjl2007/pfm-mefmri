#!/usr/bin/env python3
"""
No-MATLAB replacement for find_epi_params.m.
Creates run-level metadata sidecars used by coreg/headmotion:
  - TE.txt, TR.txt, EffectiveEchoSpacing.txt, SliceTiming.txt, PE.txt
Also updates:
  - func/qa/AllTEs.txt, func/qa/AllTRs.txt
  - func/xfms/<FuncName>/EffectiveEchoSpacing.txt (when StartSession == 1)
"""

from __future__ import annotations

import argparse
import json
import math
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
    ap.add_argument("--start-session", type=int, default=1)
    args = ap.parse_args()

    subdir = Path(args.subdir)
    func_name = args.func_name
    start_session = args.start_session
    in_root = subdir / "func" / "unprocessed" / func_name
    out_root = subdir / "func" / func_name
    qa_dir = subdir / "func" / "qa"
    xfms_dir = subdir / "func" / "xfms" / func_name

    if start_session == 1 and out_root.exists():
        shutil.rmtree(out_root)

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

            json_files = sorted(rdir.glob("Rest*.json"))
            if not json_files:
                raise FileNotFoundError(f"No Rest*.json files in {rdir}")

            te_vals = []
            last_j = None
            for jf in json_files:
                j = json.loads(jf.read_text())
                last_j = j
                te_vals.append(float(j["EchoTime"]) * 1e3)

            assert last_j is not None
            tr = float(last_j["RepetitionTime"])
            ees = float(last_j["EffectiveEchoSpacing"])
            all_ees.append(ees)

            first_j = json.loads(json_files[0].read_text())
            slice_times = [float(x) for x in first_j["SliceTiming"]]
            mb = int(first_j.get("MultibandAccelerationFactor", 1))
            slice_timing = build_slice_timing(slice_times, mb)
            phase = str(last_j.get("PhaseEncodingDirection", ""))
            pedir = phase_to_pedir(phase)

            (out_run / "TE.txt").write_text(" ".join(fmt_num(v) for v in te_vals) + "\n")
            (out_run / "TR.txt").write_text(fmt_num(tr) + "\n")
            (out_run / "EffectiveEchoSpacing.txt").write_text(fmt_num(ees) + "\n")
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
    for sdir in sessions:
        s_idx = int(sdir.name.split("_")[1])
        if s_idx < start_session:
            continue
        runs = sorted([p for p in sdir.glob("run_*") if p.is_dir()], key=lambda p: int(p.name.split("_")[1]))
        for rdir in runs:
            r_idx = int(rdir.name.split("_")[1])
            nii_files = sorted(rdir.glob("Rest_S*_R*_E*.nii.gz"))
            for e_idx, nf in enumerate(nii_files, start=1):
                size_mb = nf.stat().st_size / 1e6
                nvol = int(run(["fslnvols", str(nf)]))
                rows_fsize.append(f"{s_idx}\t{r_idx}\t{e_idx}\t{size_mb:.6f}")
                rows_nvol.append(f"{s_idx}\t{r_idx}\t{e_idx}\t{nvol}")
    fsize_tsv.write_text("session\trun\techo\tsize_mb\n" + "\n".join(rows_fsize) + ("\n" if rows_fsize else ""))
    nvol_tsv.write_text("session\trun\techo\tnvol\n" + "\n".join(rows_nvol) + ("\n" if rows_nvol else ""))


if __name__ == "__main__":
    main()

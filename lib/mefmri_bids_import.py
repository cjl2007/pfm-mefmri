#!/usr/bin/env python3
"""Convert BIDS subject inputs into RevisedMe-fMRIPipeline raw layout.

Expected output layout:
  <out_subject>/
    anat/unprocessed/T1w/T1w_*.nii.gz
    anat/unprocessed/T2w/T2w_*.nii.gz
    func/unprocessed/<func_dirname>/session_<S>/run_<R>/Rest_S<S>_R<R>_E<E>.nii.gz(+json)
    func/unprocessed/field_maps/AP_S<S>_R<R>.nii.gz(+json)
    func/unprocessed/field_maps/PA_S<S>_R<R>.nii.gz(+json)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


class ImportErrorAbort(RuntimeError):
    pass


def parse_entities(path: Path) -> Tuple[Dict[str, str], str]:
    name = path.name
    if not name.endswith(".nii.gz"):
        raise ValueError(f"Expected .nii.gz file, got: {path}")
    stem = name[:-7]
    tokens = stem.split("_")
    suffix = tokens[-1]
    entities: Dict[str, str] = {}
    for t in tokens[:-1]:
        if "-" in t:
            k, v = t.split("-", 1)
            entities[k] = v
    return entities, suffix


def subject_dir_from_bids_root(bids_root: Path, subject: str) -> Path:
    sub = subject if subject.startswith("sub-") else f"sub-{subject}"
    subdir = bids_root / sub
    if not subdir.is_dir():
        raise ImportErrorAbort(f"Missing BIDS subject dir: {subdir}")
    return subdir


def symlink_or_copy(src: Path, dst: Path, mode: str, overwrite: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if not overwrite:
            raise ImportErrorAbort(f"Refusing to overwrite existing file: {dst}")
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()

    if mode == "symlink":
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
    else:
        shutil.copy2(src, dst)


def sidecar_json_for_nii(nii: Path) -> Path:
    return nii.with_name(nii.name[:-7] + ".json")


def iter_bids_nii(subdir: Path) -> Iterable[Path]:
    for p in sorted(subdir.rglob("*.nii.gz")):
        yield p


def is_datatype_path(path: Path, dtype: str) -> bool:
    return f"/{dtype}/" in str(path)


def sort_mixed_labels(vals: Iterable[str]) -> List[str]:
    def key(v: str):
        return (0, int(v)) if v.isdigit() else (1, v)

    return sorted(set(vals), key=key)


def import_anat(
    anat_files: List[Path],
    out_subject_dir: Path,
    mode: str,
    overwrite: bool,
    warnings: List[str],
) -> Tuple[int, int]:
    t1 = []
    t2 = []
    for f in anat_files:
        entities, suffix = parse_entities(f)
        if suffix == "T1w":
            t1.append((entities, f))
        elif suffix == "T2w":
            t2.append((entities, f))

    def anat_sort_key(item):
        e, p = item
        ses = e.get("ses", "1")
        run = e.get("run", "0")
        acq = e.get("acq", "")
        return (ses, run, acq, p.name)

    t1.sort(key=anat_sort_key)
    t2.sort(key=anat_sort_key)

    t1_dir = out_subject_dir / "anat" / "unprocessed" / "T1w"
    t2_dir = out_subject_dir / "anat" / "unprocessed" / "T2w"
    t1_dir.mkdir(parents=True, exist_ok=True)
    t2_dir.mkdir(parents=True, exist_ok=True)

    for i, (_, src) in enumerate(t1, start=1):
        symlink_or_copy(src, t1_dir / f"T1w_{i}.nii.gz", mode, overwrite)
    for i, (_, src) in enumerate(t2, start=1):
        symlink_or_copy(src, t2_dir / f"T2w_{i}.nii.gz", mode, overwrite)

    if not t1:
        warnings.append("No BIDS T1w inputs found under anat/.")
    if not t2:
        warnings.append("No BIDS T2w inputs found under anat/ (pipeline can run legacy anatomical mode).")

    return len(t1), len(t2)


def run_group_key(path: Path, entities: Dict[str, str], suffix: str) -> str:
    if "run" in entities:
        return f"run-{entities['run']}"
    # If run entity is absent, derive a stable key from basename without echo.
    stem = path.name[:-7]
    stem = re.sub(r"_echo-[^_]+", "", stem)
    stem = re.sub(r"_(bold|sbref)$", "", stem)
    return stem


def map_sessions_and_runs(
    bold_files: List[Path], task: str, warnings: List[str]
) -> Tuple[Dict[str, int], Dict[Tuple[str, str], int], Dict[Tuple[str, str], List[Tuple[int, Path, Dict[str, str]]]]]:
    grouped: Dict[Tuple[str, str], List[Tuple[int, Path, Dict[str, str]]]] = defaultdict(list)
    session_labels = []

    for f in bold_files:
        entities, suffix = parse_entities(f)
        if suffix != "bold":
            continue
        if entities.get("task", "") != task:
            continue
        ses_label = entities.get("ses", "1")
        session_labels.append(ses_label)
        echo_s = entities.get("echo")
        if echo_s is None:
            warnings.append(f"Missing echo-<N> entity for bold file (skipping): {f}")
            continue
        if not echo_s.isdigit():
            warnings.append(f"Non-numeric echo entity for bold file (skipping): {f}")
            continue
        echo = int(echo_s)
        rkey = run_group_key(f, entities, suffix)
        grouped[(ses_label, rkey)].append((echo, f, entities))

    if not grouped:
        raise ImportErrorAbort(f"No task-{task} multi-echo BOLD inputs found.")

    ses_map = {ses: i for i, ses in enumerate(sort_mixed_labels(session_labels), start=1)}
    run_map: Dict[Tuple[str, str], int] = {}
    for ses in sort_mixed_labels(session_labels):
        run_keys = sorted([rk for (s, rk) in grouped.keys() if s == ses])
        for i, rk in enumerate(run_keys, start=1):
            run_map[(ses, rk)] = i

    # Per-group continuity check (warn only here; strict check happens in pipeline validator).
    for (ses, rk), entries in grouped.items():
        echoes = sorted([e for e, _, _ in entries])
        expected = list(range(1, len(echoes) + 1))
        if echoes != expected:
            warnings.append(
                f"Non-continuous echoes for BIDS group ses={ses} runkey={rk}: found {echoes}, expected {expected}"
            )

    return ses_map, run_map, grouped


def import_func(
    func_files: List[Path],
    out_subject_dir: Path,
    task: str,
    func_dirname: str,
    func_prefix: str,
    mode: str,
    overwrite: bool,
    warnings: List[str],
) -> Tuple[int, int]:
    ses_map, run_map, grouped = map_sessions_and_runs(func_files, task, warnings)
    copied_nii = 0
    copied_json = 0
    out_root = out_subject_dir / "func" / "unprocessed" / func_dirname

    for (ses_lbl, rkey), entries in sorted(grouped.items(), key=lambda x: (ses_map[x[0][0]], run_map[x[0]])):
        s = ses_map[ses_lbl]
        r = run_map[(ses_lbl, rkey)]
        run_dir = out_root / f"session_{s}" / f"run_{r}"
        run_dir.mkdir(parents=True, exist_ok=True)
        for echo, src, _ in sorted(entries, key=lambda t: t[0]):
            dst = run_dir / f"{func_prefix}_S{s}_R{r}_E{echo}.nii.gz"
            symlink_or_copy(src, dst, mode, overwrite)
            copied_nii += 1
            js = sidecar_json_for_nii(src)
            if js.exists():
                symlink_or_copy(js, run_dir / f"{func_prefix}_S{s}_R{r}_E{echo}.json", mode, overwrite)
                copied_json += 1
            else:
                warnings.append(f"Missing JSON sidecar for bold file: {src}")

    return copied_nii, copied_json


def infer_run_from_intended_for(intended_for: List[str]) -> Optional[str]:
    for rel in intended_for:
        m = re.search(r"_run-([A-Za-z0-9]+)_", rel)
        if m:
            return m.group(1)
    return None


def import_fmaps(
    fmap_files: List[Path],
    out_subject_dir: Path,
    mode: str,
    overwrite: bool,
    task: str,
    warnings: List[str],
) -> Tuple[int, int]:
    # Build minimal session/run mapping from destination func folders if they exist.
    out_func = out_subject_dir / "func" / "unprocessed"
    ses_dirs = sorted([p for p in out_func.rglob("session_*") if p.is_dir()])
    ses_map_dest: Dict[str, int] = {}
    if ses_dirs:
        # Destination sessions are already numbered; map BIDS ses labels via order fallback.
        pass

    ap_count = 0
    pa_count = 0
    out_fm = out_subject_dir / "func" / "unprocessed" / "field_maps"
    out_fm.mkdir(parents=True, exist_ok=True)

    # Group by BIDS ses/run labels for epi fieldmaps with dir entity.
    temp_groups: Dict[Tuple[str, str, str], Path] = {}
    ses_labels = set()
    run_labels_by_ses: Dict[str, set] = defaultdict(set)

    for f in fmap_files:
        entities, suffix = parse_entities(f)
        if suffix != "epi":
            continue
        dir_label = entities.get("dir", "").lower()
        if dir_label.startswith("ap"):
            pol = "AP"
        elif dir_label.startswith("pa"):
            pol = "PA"
        else:
            warnings.append(f"Skipping fmap epi without dir-AP/dir-PA entity: {f}")
            continue

        ses_lbl = entities.get("ses", "1")
        run_lbl = entities.get("run")
        if run_lbl is None:
            js = sidecar_json_for_nii(f)
            if js.exists():
                try:
                    meta = json.loads(js.read_text())
                    intended = meta.get("IntendedFor", [])
                    if isinstance(intended, list):
                        inferred = infer_run_from_intended_for(intended)
                        run_lbl = inferred
                except Exception:
                    pass
        if run_lbl is None:
            run_lbl = "1"
            warnings.append(f"No run label for fmap {f}; defaulting to run 1.")

        ses_labels.add(ses_lbl)
        run_labels_by_ses[ses_lbl].add(run_lbl)
        temp_groups[(ses_lbl, run_lbl, pol)] = f

    if not temp_groups:
        warnings.append("No BIDS fmap epi files imported (dir-AP/dir-PA not found).")
        return ap_count, pa_count

    ses_sorted = sort_mixed_labels(ses_labels)
    ses_map = {s: i for i, s in enumerate(ses_sorted, start=1)}
    run_map: Dict[Tuple[str, str], int] = {}
    for ses in ses_sorted:
        for i, run_lbl in enumerate(sort_mixed_labels(run_labels_by_ses[ses]), start=1):
            run_map[(ses, run_lbl)] = i

    for (ses_lbl, run_lbl, pol), src in sorted(
        temp_groups.items(),
        key=lambda x: (ses_map[x[0][0]], run_map[(x[0][0], x[0][1])], x[0][2]),
    ):
        s = ses_map[ses_lbl]
        r = run_map[(ses_lbl, run_lbl)]
        dst_nii = out_fm / f"{pol}_S{s}_R{r}.nii.gz"
        symlink_or_copy(src, dst_nii, mode, overwrite)
        js = sidecar_json_for_nii(src)
        if js.exists():
            symlink_or_copy(js, out_fm / f"{pol}_S{s}_R{r}.json", mode, overwrite)
        else:
            warnings.append(f"Missing JSON sidecar for fmap file: {src}")
        if pol == "AP":
            ap_count += 1
        else:
            pa_count += 1

    return ap_count, pa_count


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert BIDS subject to me-fMRI pipeline raw layout.")
    ap.add_argument("bids_root", help="Path to BIDS dataset root")
    ap.add_argument("subject", help="BIDS subject label (e.g., 06 or sub-06)")
    ap.add_argument("out_subject_dir", help="Output subject dir in pipeline layout")
    ap.add_argument("--task", default="rest", help="BIDS task label to import (default: rest)")
    ap.add_argument("--func-dirname", default="rest", help="Pipeline func subdir name (default: rest)")
    ap.add_argument("--func-prefix", default="Rest", help="Pipeline func file prefix (default: Rest)")
    ap.add_argument("--mode", choices=["symlink", "copy"], default="symlink", help="Import mode (default: symlink)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing destination files")
    args = ap.parse_args()

    bids_root = Path(args.bids_root).resolve()
    out_subject_dir = Path(args.out_subject_dir).resolve()
    subdir = subject_dir_from_bids_root(bids_root, args.subject)

    if out_subject_dir.name != (args.subject if not args.subject.startswith("sub-") else args.subject[4:]):
        # This is informational only; user may choose any output name.
        pass

    warnings: List[str] = []
    all_files = list(iter_bids_nii(subdir))
    anat_files = [p for p in all_files if is_datatype_path(p, "anat")]
    func_files = [p for p in all_files if is_datatype_path(p, "func")]
    fmap_files = [p for p in all_files if is_datatype_path(p, "fmap")]

    try:
        t1_n, t2_n = import_anat(anat_files, out_subject_dir, args.mode, args.overwrite, warnings)
        func_n, func_json_n = import_func(
            func_files,
            out_subject_dir,
            args.task,
            args.func_dirname,
            args.func_prefix,
            args.mode,
            args.overwrite,
            warnings,
        )
        ap_n, pa_n = import_fmaps(
            fmap_files,
            out_subject_dir,
            args.mode,
            args.overwrite,
            args.task,
            warnings,
        )
    except ImportErrorAbort as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print("BIDS import complete.")
    print(f"  BIDS subject dir: {subdir}")
    print(f"  Output subject dir: {out_subject_dir}")
    print(f"  Import mode: {args.mode}")
    print(f"  T1w files: {t1_n}")
    print(f"  T2w files: {t2_n}")
    print(f"  Func echoes (nii): {func_n}")
    print(f"  Func echoes (json): {func_json_n}")
    print(f"  Fieldmaps AP: {ap_n}")
    print(f"  Fieldmaps PA: {pa_n}")

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  WARNING: {w}")

    return 0


if __name__ == "__main__":
    sys.exit(main())


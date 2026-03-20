#!/usr/bin/env python3
"""Strict raw DICOM intake helper for RevisedMe-fMRIPipeline."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


class ImportAbort(RuntimeError):
    pass


@dataclass
class ProtocolConfig:
    protocol_name: str
    dcm2niix_bin: str
    func_dirname: str
    func_file_prefix: str
    expect_rest_runs_per_session: int
    expect_echoes_per_run: int
    expect_sbref_per_run: int
    expect_fmap_ap_per_session: int
    expect_fmap_pa_per_session: int
    expect_t1w_max_per_import: int
    expect_t2w_max_per_import: int
    require_t1w_if_subject_missing: bool
    t2w_optional: bool
    expect_rest_volumes: int
    expect_sbref_volumes: int
    expect_fmap_volumes: int
    expect_anat_volumes: int
    min_bytes_rest: int
    min_bytes_sbref: int
    min_bytes_fmap: int
    min_bytes_t1w: int
    min_bytes_t2w: int
    t1w_regex: str
    t2w_regex: str
    rest_regex: str
    sbref_regex: str
    fmap_ap_regex: str
    fmap_pa_regex: str
    ignore_regexes: List[str]

    @classmethod
    def from_env(cls) -> "ProtocolConfig":
        def env_int(name: str, default: int = 0) -> int:
            return int(os.environ.get(name, default))

        def env_bool(name: str, default: int = 0) -> bool:
            return bool(env_int(name, default))

        ignore_serialized = os.environ.get("IMPORT_IGNORE_REGEXES_SERIALIZED", "")
        ignore_regexes = [line for line in ignore_serialized.splitlines() if line.strip()]
        return cls(
            protocol_name=os.environ.get("IMPORT_PROTOCOL_NAME", "unknown_protocol"),
            dcm2niix_bin=os.environ.get("IMPORT_DCM2NIIX_BIN", "dcm2niix"),
            func_dirname=os.environ.get("FUNC_DIRNAME", "rest"),
            func_file_prefix=os.environ.get("FUNC_FILE_PREFIX", "Rest"),
            expect_rest_runs_per_session=env_int("IMPORT_EXPECT_REST_RUNS_PER_SESSION"),
            expect_echoes_per_run=env_int("IMPORT_EXPECT_ECHOES_PER_RUN"),
            expect_sbref_per_run=env_int("IMPORT_EXPECT_SBREF_PER_RUN"),
            expect_fmap_ap_per_session=env_int("IMPORT_EXPECT_FMAP_AP_PER_SESSION"),
            expect_fmap_pa_per_session=env_int("IMPORT_EXPECT_FMAP_PA_PER_SESSION"),
            expect_t1w_max_per_import=env_int("IMPORT_EXPECT_T1W_MAX_PER_IMPORT"),
            expect_t2w_max_per_import=env_int("IMPORT_EXPECT_T2W_MAX_PER_IMPORT"),
            require_t1w_if_subject_missing=env_bool("IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING"),
            t2w_optional=env_bool("IMPORT_T2W_OPTIONAL", 1),
            expect_rest_volumes=env_int("IMPORT_EXPECT_REST_VOLUMES"),
            expect_sbref_volumes=env_int("IMPORT_EXPECT_SBREF_VOLUMES"),
            expect_fmap_volumes=env_int("IMPORT_EXPECT_FMAP_VOLUMES"),
            expect_anat_volumes=env_int("IMPORT_EXPECT_ANAT_VOLUMES"),
            min_bytes_rest=env_int("IMPORT_MIN_BYTES_REST"),
            min_bytes_sbref=env_int("IMPORT_MIN_BYTES_SBREF"),
            min_bytes_fmap=env_int("IMPORT_MIN_BYTES_FMAP"),
            min_bytes_t1w=env_int("IMPORT_MIN_BYTES_T1W"),
            min_bytes_t2w=env_int("IMPORT_MIN_BYTES_T2W"),
            t1w_regex=os.environ.get("IMPORT_T1W_REGEX", ""),
            t2w_regex=os.environ.get("IMPORT_T2W_REGEX", ""),
            rest_regex=os.environ.get("IMPORT_REST_REGEX", ""),
            sbref_regex=os.environ.get("IMPORT_SBREF_REGEX", ""),
            fmap_ap_regex=os.environ.get("IMPORT_FMAP_AP_REGEX", ""),
            fmap_pa_regex=os.environ.get("IMPORT_FMAP_PA_REGEX", ""),
            ignore_regexes=ignore_regexes,
        )


@dataclass
class SeriesRecord:
    nii_path: Path
    json_path: Path
    basename: str
    series_description: str
    protocol_name: str
    image_type: List[str]
    series_number: int
    acquisition_time: str
    echo_number: Optional[int]
    echo_time: Optional[float]
    phase_encoding_direction: str
    n_volumes: int
    file_size: int
    filecount_images: Optional[int]
    classification: str = "unclassified"
    classification_reason: str = ""


class Reporter:
    def __init__(self, report_path: Path):
        self.report_path = report_path
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = report_path.open("w", encoding="utf-8")

    def log(self, msg: str = "") -> None:
        print(msg, flush=True)
        self._fh.write(msg + "\n")
        self._fh.flush()

    def close(self) -> None:
        self._fh.close()


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Import raw scanner-export DICOMs into me-fMRI pipeline layout.")
    ap.add_argument("raw_dicom_dir", help="Raw scanner-export DICOM directory")
    ap.add_argument("subject_dir", help="Destination subject directory")
    ap.add_argument("--session", required=True, type=int, help="Session number to assign")
    ap.add_argument("--config-file", required=True, help="Config file path used for this import")
    ap.add_argument("--dry-run", action="store_true", help="Validate and report without copying files into pipeline layout")
    return ap.parse_args()


def parse_filecount(path: Path) -> Dict[int, int]:
    mapping: Dict[int, int] = {}
    if not path.is_file():
        return mapping
    line_re = re.compile(r"^\s*([0-9]+)\s*,.*?,\s*S([0-9]+)\s*,")
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = line_re.match(line)
        if m:
            mapping[int(m.group(2))] = int(m.group(1))
    return mapping


def run_dcm2niix(raw_dir: Path, out_dir: Path, dcm2niix_bin: str, reporter: Reporter) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [dcm2niix_bin, "-z", "y", "-b", "y", "-ba", "y", "-o", str(out_dir), str(raw_dir)]
    reporter.log(f"[import] Running dcm2niix: {' '.join(cmd)}")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        reporter.log(f"[dcm2niix] {line.rstrip()}")
    rc = proc.wait()
    if rc != 0:
        raise ImportAbort(f"dcm2niix failed with exit code {rc}")


def iter_converted_series(staging_dir: Path) -> Iterable[Tuple[Path, Path]]:
    for nii_path in sorted(staging_dir.glob("*.nii.gz")):
        json_path = nii_path.with_name(nii_path.name[:-7] + ".json")
        if json_path.is_file():
            yield nii_path, json_path


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ImportAbort(f"Failed to parse JSON {path}: {exc}") from exc


def parse_nifti_volumes(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:
        header = fh.read(348)
    if len(header) < 48:
        raise ImportAbort(f"NIfTI header too short: {path}")

    sizeof_hdr = struct.unpack("<I", header[:4])[0]
    endian = "<"
    if sizeof_hdr != 348:
        sizeof_hdr_be = struct.unpack(">I", header[:4])[0]
        if sizeof_hdr_be != 348:
            raise ImportAbort(f"Invalid NIfTI header size for {path}")
        endian = ">"
    dims = struct.unpack(endian + "8h", header[40:56])
    ndim = max(dims[0], 0)
    if ndim >= 4:
        return max(int(dims[4]), 1)
    return 1


def parse_series_number(value: object, fallback_name: str) -> int:
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    m = re.search(r"_([0-9]+)_[0-9]+(?:_ph)?\.nii\.gz$", fallback_name)
    if m:
        return int(m.group(1))
    raise ImportAbort(f"Unable to parse series number for {fallback_name}")


def build_record(nii_path: Path, json_path: Path, filecount_images: Optional[int]) -> SeriesRecord:
    meta = read_json(json_path)
    image_type = [str(v) for v in meta.get("ImageType", [])]
    return SeriesRecord(
        nii_path=nii_path,
        json_path=json_path,
        basename=nii_path.name,
        series_description=str(meta.get("SeriesDescription", "")),
        protocol_name=str(meta.get("ProtocolName", "")),
        image_type=image_type,
        series_number=parse_series_number(meta.get("SeriesNumber"), nii_path.name),
        acquisition_time=str(meta.get("AcquisitionTime", "")),
        echo_number=int(meta["EchoNumber"]) if str(meta.get("EchoNumber", "")).isdigit() else None,
        echo_time=float(meta["EchoTime"]) if "EchoTime" in meta else None,
        phase_encoding_direction=str(meta.get("PhaseEncodingDirection", "")),
        n_volumes=parse_nifti_volumes(nii_path),
        file_size=nii_path.stat().st_size,
        filecount_images=filecount_images,
    )


def classify_record(record: SeriesRecord, cfg: ProtocolConfig) -> Tuple[str, str]:
    desc = record.series_description
    for pat in cfg.ignore_regexes:
        if re.search(pat, desc):
            return "ignored", f"matched ignore regex: {pat}"
    if record.basename.endswith("_ph.nii.gz") or "P" in record.image_type:
        return "ignored", "phase image"
    if re.search(cfg.t1w_regex, desc):
        return "t1w", "matched T1w regex"
    if re.search(cfg.t2w_regex, desc):
        return "t2w", "matched T2w regex"
    if re.search(cfg.fmap_ap_regex, desc):
        return "fmap_ap", "matched AP fieldmap regex"
    if re.search(cfg.fmap_pa_regex, desc):
        return "fmap_pa", "matched PA fieldmap regex"
    if re.search(cfg.sbref_regex, desc):
        return "sbref", "matched SBref regex"
    if re.search(cfg.rest_regex, desc):
        return "rest", "matched rest regex"
    return "unknown", "no classification rule matched"


def summarize_record(record: SeriesRecord) -> str:
    filecount = f" dicom_images={record.filecount_images}" if record.filecount_images is not None else ""
    return (
        f"series={record.series_number} desc='{record.series_description}' file='{record.basename}' "
        f"echo={record.echo_number} vols={record.n_volumes} bytes={record.file_size}{filecount}"
    )


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise ImportAbort(message)


def existing_anat_count(subject_dir: Path, kind: str) -> int:
    target = subject_dir / "anat" / "unprocessed" / kind
    prefix = f"{kind}_"
    if not target.is_dir():
        return 0
    return len(sorted(target.glob(f"{prefix}*.nii.gz")))


def next_anat_index(subject_dir: Path, kind: str) -> int:
    return existing_anat_count(subject_dir, kind) + 1


def validate_singletons(records: Sequence[SeriesRecord], expected_max: int, role_name: str) -> None:
    ensure(len(records) <= expected_max, f"Expected at most {expected_max} {role_name} series, found {len(records)}")
    series_numbers = {r.series_number for r in records}
    ensure(len(series_numbers) == len(records), f"Duplicate {role_name} series numbers detected")


def validate_record_guardrails(record: SeriesRecord, min_bytes: int, expect_vols: int, role_name: str) -> None:
    ensure(record.file_size >= min_bytes, f"{role_name} file too small: {record.nii_path} ({record.file_size} bytes < {min_bytes})")
    ensure(record.n_volumes == expect_vols, f"{role_name} expected {expect_vols} volume(s), found {record.n_volumes}: {record.nii_path}")


def group_multi_echo(records: Sequence[SeriesRecord], expect_echoes: int, expect_vols: int, min_bytes: int, role_name: str) -> List[List[SeriesRecord]]:
    groups: Dict[int, List[SeriesRecord]] = {}
    for rec in records:
        groups.setdefault(rec.series_number, []).append(rec)
    sorted_groups: List[List[SeriesRecord]] = []
    for series_number, group in sorted(groups.items()):
        group_sorted = sorted(group, key=lambda r: (r.echo_number or 999, r.basename))
        echoes = [r.echo_number for r in group_sorted]
        ensure(all(e is not None for e in echoes), f"{role_name} series {series_number} missing EchoNumber")
        expect = list(range(1, expect_echoes + 1))
        ensure([int(e) for e in echoes] == expect, f"{role_name} series {series_number} echo mismatch: found {echoes}, expected {expect}")
        descs = {r.series_description for r in group_sorted}
        ensure(len(descs) == 1, f"{role_name} series {series_number} has inconsistent descriptions: {sorted(descs)}")
        for rec in group_sorted:
            ensure(rec.file_size >= min_bytes, f"{role_name} file too small: {rec.nii_path} ({rec.file_size} bytes < {min_bytes})")
            ensure(rec.n_volumes == expect_vols, f"{role_name} expected {expect_vols} volume(s), found {rec.n_volumes}: {rec.nii_path}")
        sorted_groups.append(group_sorted)
    return sorted(sorted_groups, key=lambda group: (group[0].series_number, group[0].acquisition_time))


def format_group(group: Sequence[SeriesRecord]) -> str:
    first = group[0]
    echoes = ",".join(str(r.echo_number) for r in group)
    filecount = f" dicom_images={first.filecount_images}" if first.filecount_images is not None else ""
    return (
        f"series={first.series_number} desc='{first.series_description}' "
        f"echoes=[{echoes}] vols={group[0].n_volumes} acquisition={first.acquisition_time}{filecount}"
    )


def prepare_output_paths(subject_dir: Path, session: int, cfg: ProtocolConfig) -> Dict[str, Path]:
    import_dir = subject_dir / "import"
    staging_dir = import_dir / "dcm2niix"
    report_path = import_dir / "import_report.txt"
    manifest_path = import_dir / "import_manifest.tsv"
    config_snapshot = import_dir / "import_config_snapshot.sh"
    func_session_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}"
    return {
        "import_dir": import_dir,
        "staging_dir": staging_dir,
        "report_path": report_path,
        "manifest_path": manifest_path,
        "config_snapshot": config_snapshot,
        "func_session_dir": func_session_dir,
    }


def ensure_dest_absent(path: Path) -> None:
    ensure(not path.exists(), f"Refusing to overwrite existing path: {path}")


def copy_with_sidecar(src_nii: Path, src_json: Path, dst_nii: Path, dst_json: Path, dry_run: bool, reporter: Reporter) -> None:
    ensure_dest_absent(dst_nii)
    ensure_dest_absent(dst_json)
    reporter.log(f"[copy] {src_nii} -> {dst_nii}")
    reporter.log(f"[copy] {src_json} -> {dst_json}")
    if dry_run:
        return
    dst_nii.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_nii, dst_nii)
    shutil.copy2(src_json, dst_json)


def write_manifest(
    manifest_path: Path,
    rows: Sequence[Tuple[str, int, str, str, str, str, int, int, str]],
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "classification",
                "series_number",
                "series_description",
                "source_nii",
                "source_json",
                "dest_nii",
                "n_volumes",
                "file_size_bytes",
                "reason",
            ]
        )
        writer.writerows(rows)


def snapshot_config(config_file: Path, snapshot_path: Path, dry_run: bool) -> None:
    if dry_run:
        return
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_file, snapshot_path)


def main() -> int:
    args = parse_args()
    cfg = ProtocolConfig.from_env()

    raw_dir = Path(args.raw_dicom_dir).resolve()
    subject_dir = Path(args.subject_dir).resolve()
    config_file = Path(args.config_file).resolve()
    session = args.session

    ensure(raw_dir.is_dir(), f"Missing raw DICOM dir: {raw_dir}")
    ensure(session >= 1, f"Session must be >= 1, got {session}")

    paths = prepare_output_paths(subject_dir, session, cfg)
    reporter = Reporter(paths["report_path"])
    try:
        reporter.log("[import] Strict raw import started")
        reporter.log(f"[import] Raw DICOM dir: {raw_dir}")
        reporter.log(f"[import] Subject dir: {subject_dir}")
        reporter.log(f"[import] Session: {session}")
        reporter.log(f"[import] Protocol: {cfg.protocol_name}")
        reporter.log(f"[import] Config file: {config_file}")
        reporter.log(f"[import] Dry run: {int(args.dry_run)}")

        if subject_dir.exists():
            reporter.log(f"[import] Subject dir exists: {subject_dir}")
        else:
            reporter.log(f"[import] Subject dir will be created: {subject_dir}")

        ensure(not paths["func_session_dir"].exists(), f"Target session already exists: {paths['func_session_dir']}")
        if paths["staging_dir"].exists():
            raise ImportAbort(f"Import staging dir already exists: {paths['staging_dir']}")

        if not args.dry_run:
            subject_dir.mkdir(parents=True, exist_ok=True)
        run_dcm2niix(raw_dir, paths["staging_dir"], cfg.dcm2niix_bin, reporter)

        filecount_map = parse_filecount(raw_dir / "FILECOUNT2.txt")
        if filecount_map:
            reporter.log(f"[import] Parsed FILECOUNT2.txt entries: {len(filecount_map)}")
        else:
            reporter.log("[import] FILECOUNT2.txt not found or empty; continuing without DICOM image-count table")

        records: List[SeriesRecord] = []
        ignored: List[SeriesRecord] = []
        unknown: List[SeriesRecord] = []
        for nii_path, json_path in iter_converted_series(paths["staging_dir"]):
            record = build_record(nii_path, json_path, None)
            record.filecount_images = filecount_map.get(record.series_number)
            classification, reason = classify_record(record, cfg)
            record.classification = classification
            record.classification_reason = reason
            reporter.log(f"[classify] {classification.upper():9s} {summarize_record(record)} reason={reason}")
            records.append(record)
            if classification == "ignored":
                ignored.append(record)
            elif classification == "unknown":
                unknown.append(record)

        ensure(records, "No converted NIfTI series found after dcm2niix")
        if unknown:
            details = "; ".join(summarize_record(r) for r in unknown)
            raise ImportAbort(f"Unclassified converted series detected: {details}")

        t1w_records = [r for r in records if r.classification == "t1w"]
        t2w_records = [r for r in records if r.classification == "t2w"]
        fmap_ap_records = [r for r in records if r.classification == "fmap_ap"]
        fmap_pa_records = [r for r in records if r.classification == "fmap_pa"]
        sbref_records = [r for r in records if r.classification == "sbref"]
        rest_records = [r for r in records if r.classification == "rest"]

        validate_singletons(t1w_records, cfg.expect_t1w_max_per_import, "T1w")
        validate_singletons(t2w_records, cfg.expect_t2w_max_per_import, "T2w")
        for rec in t1w_records:
            validate_record_guardrails(rec, cfg.min_bytes_t1w, cfg.expect_anat_volumes, "T1w")
        for rec in t2w_records:
            validate_record_guardrails(rec, cfg.min_bytes_t2w, cfg.expect_anat_volumes, "T2w")

        if session == 1 and cfg.require_t1w_if_subject_missing:
            existing_t1w = existing_anat_count(subject_dir, "T1w")
            ensure(existing_t1w > 0 or len(t1w_records) >= 1, "Session 1 import requires at least one T1w anatomical")

        ensure(len(fmap_ap_records) == cfg.expect_fmap_ap_per_session, f"Expected {cfg.expect_fmap_ap_per_session} AP fieldmaps, found {len(fmap_ap_records)}")
        ensure(len(fmap_pa_records) == cfg.expect_fmap_pa_per_session, f"Expected {cfg.expect_fmap_pa_per_session} PA fieldmaps, found {len(fmap_pa_records)}")
        for rec in fmap_ap_records + fmap_pa_records:
            validate_record_guardrails(rec, cfg.min_bytes_fmap, cfg.expect_fmap_volumes, "Fieldmap")

        rest_groups = group_multi_echo(rest_records, cfg.expect_echoes_per_run, cfg.expect_rest_volumes, cfg.min_bytes_rest, "Rest")
        sbref_groups = group_multi_echo(sbref_records, cfg.expect_echoes_per_run, cfg.expect_sbref_volumes, cfg.min_bytes_sbref, "SBref")

        ensure(len(rest_groups) == cfg.expect_rest_runs_per_session, f"Expected {cfg.expect_rest_runs_per_session} resting-state runs, found {len(rest_groups)}")
        expected_sbref_groups = cfg.expect_rest_runs_per_session * cfg.expect_sbref_per_run
        ensure(len(sbref_groups) == expected_sbref_groups, f"Expected {expected_sbref_groups} SBref groups, found {len(sbref_groups)}")

        reporter.log(f"[validate] Rest groups: {len(rest_groups)}")
        for idx, group in enumerate(rest_groups, start=1):
            reporter.log(f"[validate]   Rest run {idx}: {format_group(group)}")
        reporter.log(f"[validate] SBref groups: {len(sbref_groups)}")
        for idx, group in enumerate(sbref_groups, start=1):
            reporter.log(f"[validate]   SBref run {idx}: {format_group(group)}")

        manifest_rows: List[Tuple[str, int, str, str, str, str, int, int, str]] = []

        for rec in ignored:
            manifest_rows.append(
                (
                    rec.classification,
                    rec.series_number,
                    rec.series_description,
                    str(rec.nii_path),
                    str(rec.json_path),
                    "",
                    rec.n_volumes,
                    rec.file_size,
                    rec.classification_reason,
                )
            )

        for rec in sorted(t1w_records, key=lambda r: (r.series_number, r.basename)):
            idx = next_anat_index(subject_dir, "T1w")
            dst_nii = subject_dir / "anat" / "unprocessed" / "T1w" / f"T1w_{idx}.nii.gz"
            dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
            manifest_rows.append(("t1w", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for rec in sorted(t2w_records, key=lambda r: (r.series_number, r.basename)):
            idx = next_anat_index(subject_dir, "T2w")
            dst_nii = subject_dir / "anat" / "unprocessed" / "T2w" / f"T2w_{idx}.nii.gz"
            dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
            manifest_rows.append(("t2w", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for idx, rec in enumerate(sorted(fmap_ap_records, key=lambda r: (r.series_number, r.acquisition_time)), start=1):
            dst_nii = subject_dir / "func" / "unprocessed" / "field_maps" / f"AP_S{session}_R{idx}.nii.gz"
            dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
            manifest_rows.append(("fmap_ap", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for idx, rec in enumerate(sorted(fmap_pa_records, key=lambda r: (r.series_number, r.acquisition_time)), start=1):
            dst_nii = subject_dir / "func" / "unprocessed" / "field_maps" / f"PA_S{session}_R{idx}.nii.gz"
            dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
            manifest_rows.append(("fmap_pa", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for run_idx, group in enumerate(rest_groups, start=1):
            run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
            for rec in group:
                echo = int(rec.echo_number or 0)
                dst_nii = run_dir / f"{cfg.func_file_prefix}_S{session}_R{run_idx}_E{echo}.nii.gz"
                dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
                manifest_rows.append(("rest", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for run_idx, group in enumerate(sbref_groups, start=1):
            run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
            for rec in group:
                echo = int(rec.echo_number or 0)
                dst_nii = run_dir / f"SBref_S{session}_R{run_idx}_E{echo}.nii.gz"
                dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter)
                manifest_rows.append(("sbref", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        write_manifest(paths["manifest_path"], manifest_rows)
        snapshot_config(config_file, paths["config_snapshot"], args.dry_run)
        reporter.log(f"[import] Wrote report: {paths['report_path']}")
        reporter.log(f"[import] Wrote manifest: {paths['manifest_path']}")
        if args.dry_run:
            reporter.log("[import] Dry run complete. No organized files were copied.")
        else:
            reporter.log("[import] Import complete.")
        return 0
    except ImportAbort as exc:
        reporter.log(f"[error] {exc}")
        return 2
    finally:
        reporter.close()


if __name__ == "__main__":
    sys.exit(main())

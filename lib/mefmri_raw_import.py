#!/usr/bin/env python3
"""Strict raw DICOM intake helper for RevisedMe-fMRIPipeline.

The importer writes functional runs and field maps into the selected task
folder under func/unprocessed/<FUNC_DIRNAME>/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shlex
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
    expect_task_runs_per_session: int
    expect_echoes_per_run: int
    expect_sbref_per_run: int
    expect_fmap_ap_per_session: int
    expect_fmap_pa_per_session: int
    expect_t1w_max_per_import: int
    expect_t2w_max_per_import: int
    require_t1w_if_subject_missing: bool
    t2w_optional: bool
    expect_task_volumes: int
    expect_sbref_volumes: int
    expect_fmap_volumes: int
    expect_anat_volumes: int
    min_bytes_task: int
    min_bytes_sbref: int
    min_bytes_fmap: int
    min_bytes_t1w: int
    min_bytes_t2w: int
    echo_dim4_policy: str
    t1w_regex: str
    t2w_regex: str
    task_regex: str
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
            expect_task_runs_per_session=env_int("IMPORT_EXPECT_TASK_RUNS_PER_SESSION", env_int("IMPORT_EXPECT_REST_RUNS_PER_SESSION")),
            expect_echoes_per_run=env_int("IMPORT_EXPECT_ECHOES_PER_RUN"),
            expect_sbref_per_run=env_int("IMPORT_EXPECT_SBREF_PER_RUN"),
            expect_fmap_ap_per_session=env_int("IMPORT_EXPECT_FMAP_AP_PER_SESSION"),
            expect_fmap_pa_per_session=env_int("IMPORT_EXPECT_FMAP_PA_PER_SESSION"),
            expect_t1w_max_per_import=env_int("IMPORT_EXPECT_T1W_MAX_PER_IMPORT"),
            expect_t2w_max_per_import=env_int("IMPORT_EXPECT_T2W_MAX_PER_IMPORT"),
            require_t1w_if_subject_missing=env_bool("IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING"),
            t2w_optional=env_bool("IMPORT_T2W_OPTIONAL", 1),
            expect_task_volumes=env_int("IMPORT_EXPECT_TASK_VOLUMES", env_int("IMPORT_EXPECT_REST_VOLUMES")),
            expect_sbref_volumes=env_int("IMPORT_EXPECT_SBREF_VOLUMES"),
            expect_fmap_volumes=env_int("IMPORT_EXPECT_FMAP_VOLUMES"),
            expect_anat_volumes=env_int("IMPORT_EXPECT_ANAT_VOLUMES"),
            min_bytes_task=env_int("IMPORT_MIN_BYTES_TASK", env_int("IMPORT_MIN_BYTES_REST")),
            min_bytes_sbref=env_int("IMPORT_MIN_BYTES_SBREF"),
            min_bytes_fmap=env_int("IMPORT_MIN_BYTES_FMAP"),
            min_bytes_t1w=env_int("IMPORT_MIN_BYTES_T1W"),
            min_bytes_t2w=env_int("IMPORT_MIN_BYTES_T2W"),
            echo_dim4_policy=os.environ.get("IMPORT_ECHO_DIM4_POLICY", "abort"),
            t1w_regex=os.environ.get("IMPORT_T1W_REGEX", ""),
            t2w_regex=os.environ.get("IMPORT_T2W_REGEX", ""),
            task_regex=os.environ.get("IMPORT_TASK_REGEX", os.environ.get("IMPORT_REST_REGEX", "")),
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
    target_n_volumes: Optional[int] = None
    noise_volumes_last: int = 0


PHASE_SUFFIX = "_phase"


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
    ap.add_argument(
        "--make-nordic-test-sessions",
        action="store_true",
        help="Deprecated pilot layout: session_1=NORDIC-denoised task data, session_2=original task data.",
    )
    ap.add_argument(
        "--nordic",
        action="store_true",
        help="Run NORDIC on detected task/rest echoes and write denoised magnitudes into the requested session; phase companions are used when available, otherwise magnitude-only NORDIC is used.",
    )
    ap.add_argument("--nordic-code-dir", default=str(Path(__file__).resolve().parent / "nordic"), help="Directory containing NIFTI_NORDIC.m")
    ap.add_argument("--nordic-matlab-bin", default="matlab", help="MATLAB executable for NORDIC")
    ap.add_argument(
        "--prepare-nordic-commands-only",
        action="store_true",
        help="Stage NORDIC outputs and write run_nordic_commands.sh, but do not launch MATLAB from Python.",
    )
    ap.add_argument(
        "--import-work-subdir",
        default="",
        help="Optional subdirectory under subject/import for staging/report/manifest outputs, useful for multi-visit imports.",
    )
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


def parse_nifti_dims(path: Path) -> Tuple[int, ...]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:
        header = fh.read(348)
    if len(header) < 56:
        raise ImportAbort(f"NIfTI header too short: {path}")

    sizeof_hdr = struct.unpack("<I", header[:4])[0]
    endian = "<"
    if sizeof_hdr != 348:
        sizeof_hdr_be = struct.unpack(">I", header[:4])[0]
        if sizeof_hdr_be != 348:
            raise ImportAbort(f"Invalid NIfTI header size for {path}")
        endian = ">"
    dims = struct.unpack(endian + "8h", header[40:56])
    ndim = max(int(dims[0]), 0)
    if ndim == 0:
        return ()
    return tuple(max(int(d), 1) for d in dims[1 : ndim + 1])


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


def is_phase_record(record: SeriesRecord) -> bool:
    return record.basename.endswith("_ph.nii.gz") or "P" in record.image_type or "PHASE" in record.image_type


def strip_phase_suffix(desc: str) -> str:
    return re.sub(r"_Pha$", "", desc)


def classify_record(record: SeriesRecord, cfg: ProtocolConfig) -> Tuple[str, str]:
    desc = record.series_description
    if is_phase_record(record):
        mag_desc = strip_phase_suffix(desc)
        if re.search(cfg.sbref_regex, mag_desc):
            return "sbref_phase", "matched SBref phase companion"
        if re.search(cfg.task_regex, mag_desc):
            return "task_phase", "matched task phase companion"
        return "ignored", "phase image"
    for pat in cfg.ignore_regexes:
        if re.search(pat, desc):
            return "ignored", f"matched ignore regex: {pat}"
    if re.search(cfg.t1w_regex, desc):
        return "t1w", "matched T1w regex"
    if re.search(cfg.t2w_regex, desc):
        return "t2w", "matched T2w regex"
    # Some scanner exports use generic or even inconsistent AP/PA descriptions.
    # For spin-echo fieldmaps, prefer the phase-encoding direction when it is available.
    if desc.startswith("SpinEchoFieldMap_"):
        ped = record.phase_encoding_direction
        if ped == "j-":
            return "fmap_ap", "generic fieldmap classified via PhaseEncodingDirection=j-"
        if ped == "j":
            return "fmap_pa", "generic fieldmap classified via PhaseEncodingDirection=j"
    if re.search(cfg.fmap_ap_regex, desc):
        return "fmap_ap", "matched AP fieldmap regex"
    if re.search(cfg.fmap_pa_regex, desc):
        return "fmap_pa", "matched PA fieldmap regex"
    if re.search(cfg.sbref_regex, desc):
        return "sbref", "matched SBref regex"
    if re.search(cfg.task_regex, desc):
        return "task", "matched task regex"
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


def prune_extra_fieldmaps(records: Sequence[SeriesRecord], expected_count: int, role_name: str, reporter: Reporter) -> List[SeriesRecord]:
    kept = list(records)
    if expected_count < 0:
        raise ImportAbort(f"{role_name} expected count must be >= 0, got {expected_count}")
    if len(kept) <= expected_count:
        return kept
    if expected_count == 0:
        details = "; ".join(summarize_record(r) for r in kept)
        raise ImportAbort(f"Unexpected {role_name} fieldmaps detected: {details}")

    if any(r.filecount_images is None for r in kept):
        details = "; ".join(summarize_record(r) for r in kept)
        raise ImportAbort(f"Expected {expected_count} {role_name} fieldmaps, found {len(kept)} and cannot disambiguate extras without FILECOUNT2.txt: {details}")

    ranked = sorted(
        kept,
        key=lambda r: (
            int(r.filecount_images or 0),
            int(r.file_size),
            r.acquisition_time,
            -int(r.series_number),
        ),
        reverse=True,
    )
    cutoff = int(ranked[expected_count - 1].filecount_images or 0)
    next_count = int(ranked[expected_count].filecount_images or 0)
    if cutoff <= next_count:
        details = "; ".join(summarize_record(r) for r in ranked)
        raise ImportAbort(f"Expected {expected_count} {role_name} fieldmaps, found {len(ranked)} with no clearly incomplete extra series: {details}")

    kept = ranked[:expected_count]
    ignored = ranked[expected_count:]
    for rec in ignored:
        reporter.log(
            f"[classify] IGNORED   {summarize_record(rec)} "
            f"reason=extra {role_name} fieldmap ignored because FILECOUNT2 indicates fewer DICOM images than retained series"
        )
        rec.classification = "ignored"
        rec.classification_reason = f"extra {role_name} fieldmap ignored due to lower FILECOUNT2 image count"
    return sorted(kept, key=lambda r: (r.series_number, r.acquisition_time, r.basename))


def group_multi_echo(
    records: Sequence[SeriesRecord],
    expect_echoes: int,
    expect_vols: int,
    min_bytes: int,
    role_name: str,
    echo_dim4_policy: str,
) -> List[List[SeriesRecord]]:
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

        if expect_vols > 0:
            for rec in group_sorted:
                if rec.n_volumes < expect_vols:
                    raise ImportAbort(f"{role_name} expected at least {expect_vols} volume(s), found {rec.n_volumes}: {rec.nii_path}")
                if rec.n_volumes > expect_vols:
                    rec.target_n_volumes = expect_vols
                    rec.noise_volumes_last = int(rec.n_volumes) - int(expect_vols)

        effective_vols = sorted({int(r.target_n_volumes or r.n_volumes) for r in group_sorted})
        if len(effective_vols) > 1:
            if echo_dim4_policy == "truncate_to_min":
                min_vol = int(min(effective_vols))
                for rec in group_sorted:
                    if rec.target_n_volumes is None or rec.target_n_volumes > min_vol:
                        rec.target_n_volumes = min_vol
            else:
                names = ", ".join(f"E{int(r.echo_number or 0)}={int(r.target_n_volumes or r.n_volumes)}" for r in group_sorted)
                raise ImportAbort(
                    f"{role_name} series {series_number} has mismatched echo dim4/timepoints ({names}). "
                    "This often indicates incomplete transfer. Re-import from source or re-run with IMPORT_ECHO_DIM4_POLICY=truncate_to_min to force truncation to min timepoints."
                )

        if expect_vols > 0 and echo_dim4_policy != "truncate_to_min":
            for rec in group_sorted:
                ensure(
                    int(rec.target_n_volumes or rec.n_volumes) == expect_vols,
                    f"{role_name} expected {expect_vols} retained volume(s), found {rec.n_volumes}: {rec.nii_path}",
                )

        sorted_groups.append(group_sorted)
    return sorted(sorted_groups, key=lambda group: (group[0].series_number, group[0].acquisition_time))


def validate_phase_groups(
    mag_groups: Sequence[Sequence[SeriesRecord]],
    phase_groups: Sequence[Sequence[SeriesRecord]],
    role_name: str,
) -> None:
    ensure(len(phase_groups) == len(mag_groups), f"Expected {len(mag_groups)} {role_name} phase groups, found {len(phase_groups)}")
    for run_idx, (mag_group, phase_group) in enumerate(zip(mag_groups, phase_groups), start=1):
        ensure(len(phase_group) == len(mag_group), f"{role_name} run {run_idx} phase echo count mismatch")
        for mag, phase in zip(mag_group, phase_group):
            echo = int(mag.echo_number or 0)
            ensure(int(phase.echo_number or 0) == echo, f"{role_name} run {run_idx} E{echo} phase echo mismatch")
            mag_retained = int(mag.target_n_volumes or mag.n_volumes)
            phase_retained = int(phase.target_n_volumes or phase.n_volumes)
            ensure(
                phase_retained == mag_retained,
                f"{role_name} run {run_idx} E{echo} phase retained-volume mismatch: mag={mag_retained} phase={phase_retained}",
            )
            if mag.echo_time is not None and phase.echo_time is not None:
                ensure(abs(float(mag.echo_time) - float(phase.echo_time)) < 1e-9, f"{role_name} run {run_idx} E{echo} phase EchoTime mismatch")
            ensure(strip_phase_suffix(phase.series_description) == mag.series_description, f"{role_name} run {run_idx} E{echo} phase description mismatch")


def format_group(group: Sequence[SeriesRecord]) -> str:
    first = group[0]
    echoes = ",".join(str(r.echo_number) for r in group)
    filecount = f" dicom_images={first.filecount_images}" if first.filecount_images is not None else ""
    return (
        f"series={first.series_number} desc='{first.series_description}' "
        f"echoes=[{echoes}] vols={group[0].n_volumes} target_vols={group[0].target_n_volumes} acquisition={first.acquisition_time}{filecount}"
    )


def prepare_output_paths(subject_dir: Path, session: int, cfg: ProtocolConfig, import_work_subdir: str = "") -> Dict[str, Path]:
    import_root = subject_dir / "import"
    import_dir = import_root / import_work_subdir if import_work_subdir else import_root
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


def copy_with_sidecar(
    src_nii: Path,
    src_json: Path,
    dst_nii: Path,
    dst_json: Path,
    dry_run: bool,
    reporter: Reporter,
    truncate_to_vols: Optional[int] = None,
) -> None:
    ensure_dest_absent(dst_nii)
    ensure_dest_absent(dst_json)
    reporter.log(f"[copy] {src_nii} -> {dst_nii}")
    reporter.log(f"[copy] {src_json} -> {dst_json}")
    if truncate_to_vols is not None:
        reporter.log(f"[copy] truncation requested: target_vols={truncate_to_vols} source={src_nii}")
    if dry_run:
        return
    dst_nii.parent.mkdir(parents=True, exist_ok=True)
    if truncate_to_vols is not None and parse_nifti_volumes(src_nii) > int(truncate_to_vols):
        extract_nifti_volumes(src_nii, dst_nii, 0, int(truncate_to_vols), "truncation")
    else:
        shutil.copy2(src_nii, dst_nii)
    shutil.copy2(src_json, dst_json)


def extract_nifti_volumes(src_nii: Path, dst_nii: Path, start: int, count: int, label: str) -> None:
    if shutil.which("fslroi") is None:
        raise ImportAbort(f"{label} requested but fslroi is not available on PATH")
    dst_nii.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["fslroi", str(src_nii), str(dst_nii), str(int(start)), str(int(count))]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise ImportAbort(f"fslroi {label} failed (exit {proc.returncode}) for {src_nii}: {proc.stdout}")


def write_json_sidecar_copy(src_json: Path, dst_json: Path, updates: dict, dry_run: bool, reporter: Reporter) -> None:
    reporter.log(f"[json] copy/update {src_json} -> {dst_json}")
    if dry_run:
        return
    data = read_json(src_json)
    data.update(updates)
    dst_json.parent.mkdir(parents=True, exist_ok=True)
    dst_json.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def save_trailing_noise_series(
    record: SeriesRecord,
    main_dst_nii: Path,
    dry_run: bool,
    reporter: Reporter,
) -> Optional[Path]:
    if int(record.noise_volumes_last) <= 0:
        return None
    start = int(record.target_n_volumes or (record.n_volumes - record.noise_volumes_last))
    count = int(record.noise_volumes_last)
    dst_dir = main_dst_nii.parent / "ImportProvenance" / "NoRFVolumes"
    dst_nii = dst_dir / f"{main_dst_nii.name[:-7]}_noRF.nii.gz"
    dst_json = json_path_for_nii(dst_nii)
    ensure_dest_absent(dst_nii)
    ensure_dest_absent(dst_json)
    reporter.log(f"[norf] save trailing no-RF/noise volumes {record.nii_path} -> {dst_nii} start={start} count={count}")
    if dry_run:
        return dst_nii
    extract_nifti_volumes(record.nii_path, dst_nii, start, count, "trailing no-RF/noise volume extraction")
    write_json_sidecar_copy(
        record.json_path,
        dst_json,
        {
            "PipelineProvenanceOnly": True,
            "PipelineNoRFVolumes": True,
            "PipelineNoRFVolumeStart": start,
            "PipelineNoRFVolumeCount": count,
            "PipelineNoRFSourceFile": str(record.nii_path),
        },
        dry_run,
        reporter,
    )
    main_dst_json = json_path_for_nii(main_dst_nii)
    if main_dst_json.is_file():
        data = read_json(main_dst_json)
        data.update(
            {
                "PipelineTrailingNoRFVolumesExcised": count,
                "PipelineTrailingNoRFVolumeStart": start,
            }
        )
        main_dst_json.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return dst_nii


def stage_nordic_source(
    record: SeriesRecord,
    stage_dir: Path,
    stem: str,
    dry_run: bool,
    reporter: Reporter,
) -> Path:
    target_vols = record.target_n_volumes
    if target_vols is None or parse_nifti_volumes(record.nii_path) <= int(target_vols):
        return record.nii_path
    staged = stage_dir / f"{stem}.nii.gz"
    ensure_dest_absent(staged)
    reporter.log(f"[nordic] stage trimmed input: {record.nii_path} -> {staged} target_vols={target_vols}")
    if not dry_run:
        extract_nifti_volumes(record.nii_path, staged, 0, int(target_vols), "NORDIC input trimming")
    return staged


def phase_path_for_magnitude(dst_nii: Path) -> Path:
    ensure(dst_nii.name.endswith(".nii.gz"), f"Expected .nii.gz path: {dst_nii}")
    return dst_nii.with_name(dst_nii.name[:-7] + f"{PHASE_SUFFIX}.nii.gz")


def json_path_for_nii(nii_path: Path) -> Path:
    ensure(nii_path.name.endswith(".nii.gz"), f"Expected .nii.gz path: {nii_path}")
    return nii_path.with_name(nii_path.name[:-7] + ".json")


def copy_phase_with_sidecar(
    phase_record: SeriesRecord,
    dst_mag_nii: Path,
    dry_run: bool,
    reporter: Reporter,
) -> Tuple[Path, Path]:
    dst_phase_nii = phase_path_for_magnitude(dst_mag_nii)
    dst_phase_json = json_path_for_nii(dst_phase_nii)
    copy_with_sidecar(phase_record.nii_path, phase_record.json_path, dst_phase_nii, dst_phase_json, dry_run, reporter, phase_record.target_n_volumes)
    return dst_phase_nii, dst_phase_json


def write_json(path: Path, data: dict, dry_run: bool, reporter: Reporter) -> None:
    reporter.log(f"[json] write {path}")
    if dry_run:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def matlab_quote(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def nordic_expected_output(out_prefix: Path) -> Path:
    return out_prefix.with_name(out_prefix.name + ".nii.gz")


def build_nordic_expr(repo_lib_dir: Path, nordic_code_dir: Path, mag: Path, phase: Optional[Path], out_prefix: Path, log_path: Path) -> str:
    phase_arg = "" if phase is None else str(phase)
    return (
        f"addpath({matlab_quote(repo_lib_dir)}); "
        f"run_nordic_mefmri({matlab_quote(mag)}, {matlab_quote(phase_arg)}, "
        f"{matlab_quote(out_prefix)}, {matlab_quote(nordic_code_dir)}, {matlab_quote(log_path)})"
    )


def build_nordic_shell_line(matlab_bin: str, repo_lib_dir: Path, nordic_code_dir: Path, mag: Path, phase: Optional[Path], out_prefix: Path, log_path: Path) -> str:
    expr = build_nordic_expr(repo_lib_dir, nordic_code_dir, mag, phase, out_prefix, log_path)
    return f"{shlex.quote(matlab_bin)} -batch {shlex.quote(expr)}"


def build_nordic_command(matlab_bin: str, repo_lib_dir: Path, nordic_code_dir: Path, mag: Path, phase: Optional[Path], out_prefix: Path, log_path: Path) -> List[str]:
    expr = build_nordic_expr(repo_lib_dir, nordic_code_dir, mag, phase, out_prefix, log_path)
    # MATLAB activation/licensing can fail when launched directly by Python's
    # subprocess on this workstation, while the same -batch call succeeds from
    # an interactive shell. Keep the shell boundary explicit and quoted.
    shell_cmd = f"exec {shlex.quote(matlab_bin)} -batch {shlex.quote(expr)}"
    return ["bash", "-lc", shell_cmd]


def run_nordic_pair(
    matlab_bin: str,
    repo_lib_dir: Path,
    nordic_code_dir: Path,
    mag: Path,
    phase: Optional[Path],
    out_prefix: Path,
    log_path: Path,
    dry_run: bool,
    reporter: Reporter,
) -> Path:
    cmd = build_nordic_command(matlab_bin, repo_lib_dir, nordic_code_dir, mag, phase, out_prefix, log_path)
    reporter.log(f"[nordic] MATLAB command: {' '.join(cmd)}")
    if dry_run:
        return nordic_expected_output(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.stdout:
        reporter.log(proc.stdout.rstrip())
    if proc.returncode != 0:
        raise ImportAbort(f"NORDIC failed with exit code {proc.returncode}: {mag}")
    out_nii = nordic_expected_output(out_prefix)
    ensure(out_nii.is_file(), f"NORDIC output missing after MATLAB completed: {out_nii}")
    ensure(parse_nifti_dims(out_nii) == parse_nifti_dims(mag), f"NORDIC output dimensions differ from input: {out_nii}")
    ensure(parse_nifti_volumes(out_nii) == parse_nifti_volumes(mag), f"NORDIC output volume count differs from input: {out_nii}")
    return out_nii


def write_nordic_json(
    src_json: Path,
    dst_json: Path,
    src_mag: Path,
    src_phase: Optional[Path],
    dry_run: bool,
    reporter: Reporter,
    *,
    source_noise_volumes_last: int = 0,
) -> None:
    data = read_json(src_json)
    updates = {
        "Denoising": "NORDIC",
        "NORDICNoiseVolumeLast": 0,
        "NORDICUsedMagnitude": True,
        "NORDICUsedPhase": src_phase is not None,
        "NORDICSourceMagnitude": str(src_mag),
        "NORDICTrailingNoRFVolumesPresentInSource": int(source_noise_volumes_last),
        "NORDICTrailingNoRFVolumesExcised": int(source_noise_volumes_last),
        "NORDICNoRFVolumesUsed": False,
    }
    if src_phase is not None:
        updates["NORDICSourcePhase"] = str(src_phase)
    data.update(updates)
    write_json(dst_json, data, dry_run, reporter)


def write_nordic_marker(
    run_dir: Path,
    rows: Sequence[Tuple[int, Path, Optional[Path], Path]],
    dry_run: bool,
    reporter: Reporter,
    *,
    used_phase: Optional[bool] = None,
    noise_volume_last: int = 0,
    excised_noise_volume_last: int = 0,
) -> None:
    if not rows:
        return
    marker = run_dir / "NORDIC_DENOISING.txt"
    if used_phase is None:
        used_phase = any(src_phase is not None for _, _, src_phase, _ in rows)
    mode = "magnitude+phase" if used_phase else "magnitude-only"
    reporter.log(f"[nordic] write marker: {marker}")
    if dry_run:
        return
    run_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "NORDIC denoising: yes",
        f"Mode: {mode}",
        "Used magnitude: yes",
        f"Used phase: {'yes' if used_phase else 'no'}",
        f"Trailing no-RF/noise volumes present in source: {'yes' if int(excised_noise_volume_last) > 0 else 'no'}",
        f"Trailing no-RF/noise volumes excised from main run: {int(excised_noise_volume_last)}",
        f"No-RF/noise volumes used by NORDIC: {'yes' if int(noise_volume_last) > 0 else 'no'}",
        f"Trailing no-RF/noise volumes used by NORDIC: {int(noise_volume_last)}",
        "",
        "Echo outputs:",
    ]
    for echo, src_mag, src_phase, out_mag in rows:
        lines.append(f"  E{echo}: {out_mag}")
        lines.append(f"    source magnitude: {src_mag}")
        if src_phase is not None:
            lines.append(f"    source phase: {src_phase}")
    marker.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_nordic_command_script(
    script_path: Path,
    qc_path: Path,
    rows: Sequence[Tuple[int, int, Path, Optional[Path], Path, str]],
    dry_run: bool,
    reporter: Reporter,
) -> None:
    reporter.log(f"[nordic] write command script: {script_path}")
    if dry_run:
        return

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"QC_PATH={shlex.quote(str(qc_path))}",
        "printf 'run\\techo\\tinput_dims\\toutput_dims\\tvolumes\\tstatus\\n' > \"$QC_PATH\"",
        "",
    ]

    for run_idx, echo, mag, phase, out_nii, command_line in rows:
        log_label = f"S1_R{run_idx}_E{echo}"
        phase_label = "<none; magnitude-only>" if phase is None else str(phase)
        lines.extend(
            [
                f"echo '[nordic-script] {log_label}'",
                f"echo '[nordic-script] magnitude: {mag}'",
                f"echo '[nordic-script] phase: {phase_label}'",
                command_line,
                (
                    "python3 - "
                    f"{shlex.quote(str(mag))} {shlex.quote(str(out_nii))} "
                    f"\"$QC_PATH\" {run_idx} {echo} <<'PY'"
                ),
                "import gzip, struct, sys",
                "mag, out, qc, run_idx, echo = sys.argv[1:6]",
                "def dims(path):",
                "    opener = gzip.open if path.endswith('.gz') else open",
                "    with opener(path, 'rb') as f:",
                "        hdr = f.read(348)",
                "    endian = '<' if struct.unpack('<I', hdr[:4])[0] == 348 else '>'",
                "    d = struct.unpack(endian + '8h', hdr[40:56])",
                "    return tuple(int(x) for x in d[1:int(d[0]) + 1])",
                "mi = dims(mag)",
                "mo = dims(out)",
                "if mi != mo:",
                "    raise SystemExit(f'NORDIC dimension mismatch: {mag} {mi} != {out} {mo}')",
                "vols = mi[3] if len(mi) > 3 else 1",
                "with open(qc, 'a', encoding='utf-8') as f:",
                "    f.write(f'{run_idx}\\t{echo}\\t' + 'x'.join(map(str, mi)) + '\\t' + 'x'.join(map(str, mo)) + f'\\t{vols}\\tOK\\n')",
                "PY",
                "",
            ]
        )

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(lines), encoding="utf-8")
    script_path.chmod(0o755)


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
    use_nordic = bool(args.nordic or args.make_nordic_test_sessions or args.prepare_nordic_commands_only)

    raw_dir = Path(args.raw_dicom_dir).resolve()
    subject_dir = Path(args.subject_dir).resolve()
    config_file = Path(args.config_file).resolve()
    session = args.session

    ensure(raw_dir.is_dir(), f"Missing raw DICOM dir: {raw_dir}")
    ensure(session >= 1, f"Session must be >= 1, got {session}")
    ensure(not args.import_work_subdir or not Path(args.import_work_subdir).is_absolute(), "--import-work-subdir must be relative")

    paths = prepare_output_paths(subject_dir, session, cfg, args.import_work_subdir)
    reporter = Reporter(paths["report_path"])
    try:
        reporter.log("[import] Strict raw import started")
        reporter.log(f"[import] Raw DICOM dir: {raw_dir}")
        reporter.log(f"[import] Subject dir: {subject_dir}")
        reporter.log(f"[import] Session: {session}")
        reporter.log(f"[import] Protocol: {cfg.protocol_name}")
        reporter.log(f"[import] Config file: {config_file}")
        reporter.log(f"[import] Dry run: {int(args.dry_run)}")
        reporter.log(f"[import] Echo dim4 policy: {cfg.echo_dim4_policy}")
        reporter.log(f"[import] NORDIC enabled: {int(use_nordic)}")
        reporter.log(f"[import] NORDIC test sessions: {int(args.make_nordic_test_sessions)}")
        if args.import_work_subdir:
            reporter.log(f"[import] Import work subdir: {args.import_work_subdir}")
        if use_nordic and not args.make_nordic_test_sessions:
            reporter.log("[import] NORDIC layout: denoised magnitudes written into requested session")
        if args.make_nordic_test_sessions:
            reporter.log("[import] NORDIC test layout: session_1=NORDIC, session_2=original")
        if use_nordic:
            reporter.log(f"[import] NORDIC code dir: {args.nordic_code_dir}")
            reporter.log(f"[import] NORDIC MATLAB bin: {args.nordic_matlab_bin}")

        if subject_dir.exists():
            reporter.log(f"[import] Subject dir exists: {subject_dir}")
        else:
            reporter.log(f"[import] Subject dir will be created: {subject_dir}")

        if args.make_nordic_test_sessions:
            for ses in (1, 2):
                target_session = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{ses}"
                ensure(not target_session.exists(), f"Target session already exists: {target_session}")
        else:
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
        task_records = [r for r in records if r.classification == "task"]
        task_phase_records = [r for r in records if r.classification == "task_phase"]
        sbref_phase_records = [r for r in records if r.classification == "sbref_phase"]

        validate_singletons(t1w_records, cfg.expect_t1w_max_per_import, "T1w")
        validate_singletons(t2w_records, cfg.expect_t2w_max_per_import, "T2w")
        for rec in t1w_records:
            validate_record_guardrails(rec, cfg.min_bytes_t1w, cfg.expect_anat_volumes, "T1w")
        for rec in t2w_records:
            validate_record_guardrails(rec, cfg.min_bytes_t2w, cfg.expect_anat_volumes, "T2w")

        if session == 1 and cfg.require_t1w_if_subject_missing:
            existing_t1w = existing_anat_count(subject_dir, "T1w")
            ensure(existing_t1w > 0 or len(t1w_records) >= 1, "Session 1 import requires at least one T1w anatomical")

        fmap_ap_records = prune_extra_fieldmaps(fmap_ap_records, cfg.expect_fmap_ap_per_session, "AP", reporter)
        fmap_pa_records = prune_extra_fieldmaps(fmap_pa_records, cfg.expect_fmap_pa_per_session, "PA", reporter)
        ensure(len(fmap_ap_records) == cfg.expect_fmap_ap_per_session, f"Expected {cfg.expect_fmap_ap_per_session} AP fieldmaps, found {len(fmap_ap_records)}")
        ensure(len(fmap_pa_records) == cfg.expect_fmap_pa_per_session, f"Expected {cfg.expect_fmap_pa_per_session} PA fieldmaps, found {len(fmap_pa_records)}")
        for rec in fmap_ap_records + fmap_pa_records:
            validate_record_guardrails(rec, cfg.min_bytes_fmap, cfg.expect_fmap_volumes, "Fieldmap")

        ensure(cfg.echo_dim4_policy in {"abort", "truncate_to_min"}, f"IMPORT_ECHO_DIM4_POLICY must be abort|truncate_to_min, got: {cfg.echo_dim4_policy}")

        task_groups = group_multi_echo(task_records, cfg.expect_echoes_per_run, cfg.expect_task_volumes, cfg.min_bytes_task, cfg.func_file_prefix, cfg.echo_dim4_policy)
        sbref_groups = group_multi_echo(sbref_records, cfg.expect_echoes_per_run, cfg.expect_sbref_volumes, cfg.min_bytes_sbref, "SBref", cfg.echo_dim4_policy)
        task_phase_groups = group_multi_echo(task_phase_records, cfg.expect_echoes_per_run, cfg.expect_task_volumes, 0, f"{cfg.func_file_prefix} phase", cfg.echo_dim4_policy)
        sbref_phase_groups = group_multi_echo(sbref_phase_records, cfg.expect_echoes_per_run, cfg.expect_sbref_volumes, 0, "SBref phase", cfg.echo_dim4_policy)

        ensure(len(task_groups) == cfg.expect_task_runs_per_session, f"Expected {cfg.expect_task_runs_per_session} {cfg.func_dirname} runs, found {len(task_groups)}")
        expected_sbref_groups = cfg.expect_task_runs_per_session * cfg.expect_sbref_per_run
        ensure(len(sbref_groups) == expected_sbref_groups, f"Expected {expected_sbref_groups} SBref groups, found {len(sbref_groups)}")
        if task_phase_groups:
            validate_phase_groups(task_groups, task_phase_groups, cfg.func_file_prefix)
        if sbref_phase_groups:
            validate_phase_groups(sbref_groups, sbref_phase_groups, "SBref")

        reporter.log(f"[validate] {cfg.func_dirname} groups: {len(task_groups)}")
        for idx, group in enumerate(task_groups, start=1):
            reporter.log(f"[validate]   {cfg.func_dirname} run {idx}: {format_group(group)}")
        reporter.log(f"[validate] {cfg.func_dirname} phase groups: {len(task_phase_groups)}")
        for idx, group in enumerate(task_phase_groups, start=1):
            reporter.log(f"[validate]   {cfg.func_dirname} phase run {idx}: {format_group(group)}")
        reporter.log(f"[validate] SBref groups: {len(sbref_groups)}")
        for idx, group in enumerate(sbref_groups, start=1):
            reporter.log(f"[validate]   SBref run {idx}: {format_group(group)}")
        reporter.log(f"[validate] SBref phase groups: {len(sbref_phase_groups)}")
        for idx, group in enumerate(sbref_phase_groups, start=1):
            reporter.log(f"[validate]   SBref phase run {idx}: {format_group(group)}")

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
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
            manifest_rows.append(("t1w", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        for rec in sorted(t2w_records, key=lambda r: (r.series_number, r.basename)):
            idx = next_anat_index(subject_dir, "T2w")
            dst_nii = subject_dir / "anat" / "unprocessed" / "T2w" / f"T2w_{idx}.nii.gz"
            dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
            copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
            manifest_rows.append(("t2w", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        fmap_sessions = (1, 2) if args.make_nordic_test_sessions else (session,)
        for fmap_session in fmap_sessions:
            for idx, rec in enumerate(sorted(fmap_ap_records, key=lambda r: (r.series_number, r.acquisition_time)), start=1):
                dst_nii = subject_dir / "func" / "unprocessed" / cfg.func_dirname / "field_maps" / f"AP_S{fmap_session}_R{idx}.nii.gz"
                dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                manifest_rows.append(("fmap_ap", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

            for idx, rec in enumerate(sorted(fmap_pa_records, key=lambda r: (r.series_number, r.acquisition_time)), start=1):
                dst_nii = subject_dir / "func" / "unprocessed" / cfg.func_dirname / "field_maps" / f"PA_S{fmap_session}_R{idx}.nii.gz"
                dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                manifest_rows.append(("fmap_pa", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))

        if args.make_nordic_test_sessions:
            original_session = 2
            nordic_session = 1
            repo_lib_dir = Path(__file__).resolve().parent
            nordic_code_dir = Path(args.nordic_code_dir).resolve()
            ensure(nordic_code_dir.is_dir(), f"Missing NORDIC code dir: {nordic_code_dir}")
            ensure((nordic_code_dir / "NIFTI_NORDIC.m").is_file(), f"Missing NIFTI_NORDIC.m in {nordic_code_dir}")
            nordic_log_dir = subject_dir / "import" / "nordic_logs"
            nordic_qc_rows = [
                "run\techo\toriginal_magnitude\tphase\tnordic_magnitude\tinput_dims\toutput_dims\tvolumes"
            ]
            nordic_command_rows: List[Tuple[int, int, Path, Optional[Path], Path, str]] = []

            for run_idx, mag_group in enumerate(task_groups, start=1):
                phase_group = task_phase_groups[run_idx - 1] if run_idx - 1 < len(task_phase_groups) else []
                ensure(not phase_group or len(phase_group) == len(mag_group), f"{cfg.func_file_prefix} run {run_idx} phase echo count mismatch")
                original_run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{original_session}" / f"run_{run_idx}"
                nordic_run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{nordic_session}" / f"run_{run_idx}"
                nordic_marker_rows: List[Tuple[int, Path, Optional[Path], Path]] = []
                max_noise_last = max((int(rec.noise_volumes_last) for rec in mag_group), default=0)
                for echo_idx, rec in enumerate(mag_group):
                    phase_rec = phase_group[echo_idx] if phase_group else None
                    echo = int(rec.echo_number or 0)
                    original_nii = original_run_dir / f"{cfg.func_file_prefix}_S{original_session}_R{run_idx}_E{echo}.nii.gz"
                    original_json = json_path_for_nii(original_nii)
                    copy_with_sidecar(rec.nii_path, rec.json_path, original_nii, original_json, args.dry_run, reporter, rec.target_n_volumes)
                    no_rf_nii = save_trailing_noise_series(rec, original_nii, args.dry_run, reporter)
                    if no_rf_nii is not None:
                        manifest_rows.append(("task_noRF", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(no_rf_nii), rec.noise_volumes_last, rec.file_size, "trailing no-RF/noise volumes saved for provenance only"))
                    original_phase_nii = None
                    if phase_rec is not None:
                        original_phase_nii, _ = copy_phase_with_sidecar(phase_rec, original_nii, args.dry_run, reporter)
                        phase_no_rf_nii = save_trailing_noise_series(phase_rec, original_phase_nii, args.dry_run, reporter)
                        if phase_no_rf_nii is not None:
                            manifest_rows.append(("task_phase_noRF", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(phase_no_rf_nii), phase_rec.noise_volumes_last, phase_rec.file_size, "trailing no-RF/noise phase volumes saved for provenance only"))
                    manifest_rows.append(("task_original", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(original_nii), rec.n_volumes, rec.file_size, rec.classification_reason))
                    if phase_rec is not None and original_phase_nii is not None:
                        manifest_rows.append(("task_phase_original", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(original_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))

                    nordic_nii = nordic_run_dir / f"{cfg.func_file_prefix}_S{nordic_session}_R{run_idx}_E{echo}.nii.gz"
                    nordic_json = json_path_for_nii(nordic_nii)
                    nordic_phase_nii = None
                    if phase_rec is not None:
                        nordic_phase_nii, _ = copy_phase_with_sidecar(phase_rec, nordic_nii, args.dry_run, reporter)
                    out_prefix = nordic_nii.with_name(nordic_nii.name[:-7])
                    nordic_log = nordic_log_dir / f"NORDIC_S{nordic_session}_R{run_idx}_E{echo}.log"
                    if args.prepare_nordic_commands_only:
                        command_line = build_nordic_shell_line(
                            args.nordic_matlab_bin,
                            repo_lib_dir,
                            nordic_code_dir,
                            original_nii,
                            original_phase_nii,
                            out_prefix,
                            nordic_log,
                        )
                        reporter.log(f"[nordic] queued shell command: {command_line}")
                        nordic_command_rows.append((run_idx, echo, original_nii, original_phase_nii, nordic_nii, command_line))
                        out_nii = nordic_expected_output(out_prefix)
                    else:
                        out_nii = run_nordic_pair(
                            args.nordic_matlab_bin,
                            repo_lib_dir,
                            nordic_code_dir,
                            original_nii,
                            original_phase_nii,
                            out_prefix,
                            nordic_log,
                            args.dry_run,
                            reporter,
                        )
                        if not args.dry_run:
                            ensure(out_nii == nordic_nii, f"Unexpected NORDIC output path: {out_nii} != {nordic_nii}")
                    write_nordic_json(rec.json_path, nordic_json, original_nii, original_phase_nii, args.dry_run, reporter, source_noise_volumes_last=rec.noise_volumes_last)
                    manifest_rows.append(("task_nordic", rec.series_number, rec.series_description, str(original_nii), str(original_json), str(nordic_nii), rec.n_volumes, rec.file_size, "NORDIC denoised magnitude"))
                    if phase_rec is not None and nordic_phase_nii is not None:
                        manifest_rows.append(("task_phase_nordic", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(nordic_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))
                    nordic_qc_rows.append(
                        f"{run_idx}\t{echo}\t{original_nii}\t{original_phase_nii if original_phase_nii is not None else 'MAG_ONLY'}\t{nordic_nii}\t"
                        f"{'x'.join(map(str, parse_nifti_dims(rec.nii_path)))}\t"
                        f"{'PENDING' if args.prepare_nordic_commands_only else ('DRYRUN' if args.dry_run else 'x'.join(map(str, parse_nifti_dims(nordic_nii))))}\t{rec.n_volumes}"
                    )
                    nordic_marker_rows.append((echo, original_nii, original_phase_nii, nordic_nii))
                write_nordic_marker(nordic_run_dir, nordic_marker_rows, args.dry_run, reporter, excised_noise_volume_last=max_noise_last)

            for branch_session in (nordic_session, original_session):
                for run_idx, group in enumerate(sbref_groups, start=1):
                    phase_group = sbref_phase_groups[run_idx - 1] if sbref_phase_groups else []
                    run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{branch_session}" / f"run_{run_idx}"
                    for echo_idx, rec in enumerate(group):
                        echo = int(rec.echo_number or 0)
                        dst_nii = run_dir / f"SBref_S{branch_session}_R{run_idx}_E{echo}.nii.gz"
                        dst_json = json_path_for_nii(dst_nii)
                        copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                        manifest_rows.append(("sbref", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))
                        if phase_group:
                            phase_rec = phase_group[echo_idx]
                            dst_phase_nii, _ = copy_phase_with_sidecar(phase_rec, dst_nii, args.dry_run, reporter)
                            manifest_rows.append(("sbref_phase", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(dst_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))

            qc_path = subject_dir / "import" / "nordic_import_qc.tsv"
            if args.prepare_nordic_commands_only:
                write_nordic_command_script(
                    subject_dir / "import" / "run_nordic_commands.sh",
                    subject_dir / "import" / "nordic_execution_qc.tsv",
                    nordic_command_rows,
                    args.dry_run,
                    reporter,
                )
            reporter.log(f"[nordic] write QC: {qc_path}")
            if not args.dry_run:
                qc_path.parent.mkdir(parents=True, exist_ok=True)
                qc_path.write_text("\n".join(nordic_qc_rows) + "\n", encoding="utf-8")
        elif use_nordic:
            repo_lib_dir = Path(__file__).resolve().parent
            nordic_code_dir = Path(args.nordic_code_dir).resolve()
            ensure(nordic_code_dir.is_dir(), f"Missing NORDIC code dir: {nordic_code_dir}")
            ensure((nordic_code_dir / "NIFTI_NORDIC.m").is_file(), f"Missing NIFTI_NORDIC.m in {nordic_code_dir}")
            nordic_log_dir = paths["import_dir"] / "nordic_logs"
            nordic_qc_rows = [
                "session\trun\techo\tsource_magnitude\tphase\tnordic_magnitude\tinput_dims\toutput_dims\tvolumes"
            ]
            nordic_command_rows: List[Tuple[int, int, Path, Optional[Path], Path, str]] = []

            for run_idx, mag_group in enumerate(task_groups, start=1):
                phase_group = task_phase_groups[run_idx - 1] if run_idx - 1 < len(task_phase_groups) else []
                ensure(not phase_group or len(phase_group) == len(mag_group), f"{cfg.func_file_prefix} run {run_idx} phase echo count mismatch")
                run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
                stage_dir = paths["import_dir"] / "nordic_inputs" / f"session_{session}" / f"run_{run_idx}"
                nordic_marker_rows: List[Tuple[int, Path, Optional[Path], Path]] = []
                max_noise_last = max((int(rec.noise_volumes_last) for rec in mag_group), default=0)
                for echo_idx, rec in enumerate(mag_group):
                    phase_rec = phase_group[echo_idx] if phase_group else None
                    echo = int(rec.echo_number or 0)
                    nordic_nii = run_dir / f"{cfg.func_file_prefix}_S{session}_R{run_idx}_E{echo}.nii.gz"
                    nordic_json = json_path_for_nii(nordic_nii)
                    no_rf_nii = save_trailing_noise_series(rec, nordic_nii, args.dry_run, reporter)
                    if no_rf_nii is not None:
                        manifest_rows.append(("task_noRF", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(no_rf_nii), rec.noise_volumes_last, rec.file_size, "trailing no-RF/noise volumes saved for provenance only"))
                    nordic_phase_nii = None
                    if phase_rec is not None:
                        nordic_phase_nii, _ = copy_phase_with_sidecar(phase_rec, nordic_nii, args.dry_run, reporter)
                        phase_no_rf_nii = save_trailing_noise_series(phase_rec, nordic_phase_nii, args.dry_run, reporter)
                        if phase_no_rf_nii is not None:
                            manifest_rows.append(("task_phase_noRF", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(phase_no_rf_nii), phase_rec.noise_volumes_last, phase_rec.file_size, "trailing no-RF/noise phase volumes saved for provenance only"))
                    nordic_input_nii = stage_nordic_source(rec, stage_dir, f"{cfg.func_file_prefix}_S{session}_R{run_idx}_E{echo}", args.dry_run, reporter)
                    nordic_phase_input_nii = None
                    if phase_rec is not None:
                        nordic_phase_input_nii = stage_nordic_source(phase_rec, stage_dir, f"{cfg.func_file_prefix}_S{session}_R{run_idx}_E{echo}_phase", args.dry_run, reporter)
                    out_prefix = nordic_nii.with_name(nordic_nii.name[:-7])
                    nordic_log = nordic_log_dir / f"NORDIC_S{session}_R{run_idx}_E{echo}.log"
                    if args.prepare_nordic_commands_only:
                        command_line = build_nordic_shell_line(
                            args.nordic_matlab_bin,
                            repo_lib_dir,
                            nordic_code_dir,
                            nordic_input_nii,
                            nordic_phase_input_nii,
                            out_prefix,
                            nordic_log,
                        )
                        reporter.log(f"[nordic] queued shell command: {command_line}")
                        nordic_command_rows.append((run_idx, echo, nordic_input_nii, nordic_phase_input_nii, nordic_nii, command_line))
                        out_nii = nordic_expected_output(out_prefix)
                    else:
                        out_nii = run_nordic_pair(
                            args.nordic_matlab_bin,
                            repo_lib_dir,
                            nordic_code_dir,
                            nordic_input_nii,
                            nordic_phase_input_nii,
                            out_prefix,
                            nordic_log,
                            args.dry_run,
                            reporter,
                        )
                        if not args.dry_run:
                            ensure(out_nii == nordic_nii, f"Unexpected NORDIC output path: {out_nii} != {nordic_nii}")
                    write_nordic_json(rec.json_path, nordic_json, nordic_input_nii, nordic_phase_input_nii, args.dry_run, reporter, source_noise_volumes_last=rec.noise_volumes_last)
                    manifest_rows.append(("task_nordic", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(nordic_nii), rec.n_volumes, rec.file_size, "NORDIC denoised magnitude"))
                    if phase_rec is not None and nordic_phase_nii is not None:
                        manifest_rows.append(("task_phase", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(nordic_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))
                    nordic_qc_rows.append(
                        f"{session}\t{run_idx}\t{echo}\t{nordic_input_nii}\t{nordic_phase_input_nii if nordic_phase_input_nii is not None else 'MAG_ONLY'}\t{nordic_nii}\t"
                        f"{'x'.join(map(str, parse_nifti_dims(rec.nii_path)))}\t"
                        f"{'PENDING' if args.prepare_nordic_commands_only else ('DRYRUN' if args.dry_run else 'x'.join(map(str, parse_nifti_dims(nordic_nii))))}\t{rec.n_volumes}"
                    )
                    nordic_marker_rows.append((echo, nordic_input_nii, nordic_phase_input_nii, nordic_nii))
                write_nordic_marker(run_dir, nordic_marker_rows, args.dry_run, reporter, excised_noise_volume_last=max_noise_last)

            for run_idx, group in enumerate(sbref_groups, start=1):
                phase_group = sbref_phase_groups[run_idx - 1] if sbref_phase_groups else []
                run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
                for echo_idx, rec in enumerate(group):
                    echo = int(rec.echo_number or 0)
                    dst_nii = run_dir / f"SBref_S{session}_R{run_idx}_E{echo}.nii.gz"
                    dst_json = json_path_for_nii(dst_nii)
                    copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                    manifest_rows.append(("sbref", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))
                    if phase_group:
                        phase_rec = phase_group[echo_idx]
                        dst_phase_nii, _ = copy_phase_with_sidecar(phase_rec, dst_nii, args.dry_run, reporter)
                        manifest_rows.append(("sbref_phase", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(dst_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))

            qc_path = paths["import_dir"] / "nordic_import_qc.tsv"
            if args.prepare_nordic_commands_only:
                write_nordic_command_script(
                    paths["import_dir"] / "run_nordic_commands.sh",
                    paths["import_dir"] / "nordic_execution_qc.tsv",
                    nordic_command_rows,
                    args.dry_run,
                    reporter,
                )
            reporter.log(f"[nordic] write QC: {qc_path}")
            if not args.dry_run:
                qc_path.parent.mkdir(parents=True, exist_ok=True)
                qc_path.write_text("\n".join(nordic_qc_rows) + "\n", encoding="utf-8")
        else:
            for run_idx, group in enumerate(task_groups, start=1):
                phase_group = task_phase_groups[run_idx - 1] if task_phase_groups else []
                run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
                for echo_idx, rec in enumerate(group):
                    echo = int(rec.echo_number or 0)
                    dst_nii = run_dir / f"{cfg.func_file_prefix}_S{session}_R{run_idx}_E{echo}.nii.gz"
                    dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                    copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                    no_rf_nii = save_trailing_noise_series(rec, dst_nii, args.dry_run, reporter)
                    if no_rf_nii is not None:
                        manifest_rows.append(("task_noRF", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(no_rf_nii), rec.noise_volumes_last, rec.file_size, "trailing no-RF/noise volumes saved for provenance only"))
                    manifest_rows.append(("task", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))
                    if phase_group:
                        phase_rec = phase_group[echo_idx]
                        dst_phase_nii, _ = copy_phase_with_sidecar(phase_rec, dst_nii, args.dry_run, reporter)
                        phase_no_rf_nii = save_trailing_noise_series(phase_rec, dst_phase_nii, args.dry_run, reporter)
                        if phase_no_rf_nii is not None:
                            manifest_rows.append(("task_phase_noRF", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(phase_no_rf_nii), phase_rec.noise_volumes_last, phase_rec.file_size, "trailing no-RF/noise phase volumes saved for provenance only"))
                        manifest_rows.append(("task_phase", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(dst_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))

            for run_idx, group in enumerate(sbref_groups, start=1):
                phase_group = sbref_phase_groups[run_idx - 1] if sbref_phase_groups else []
                run_dir = subject_dir / "func" / "unprocessed" / cfg.func_dirname / f"session_{session}" / f"run_{run_idx}"
                for echo_idx, rec in enumerate(group):
                    echo = int(rec.echo_number or 0)
                    dst_nii = run_dir / f"SBref_S{session}_R{run_idx}_E{echo}.nii.gz"
                    dst_json = dst_nii.with_name(dst_nii.name[:-7] + ".json")
                    copy_with_sidecar(rec.nii_path, rec.json_path, dst_nii, dst_json, args.dry_run, reporter, rec.target_n_volumes)
                    manifest_rows.append(("sbref", rec.series_number, rec.series_description, str(rec.nii_path), str(rec.json_path), str(dst_nii), rec.n_volumes, rec.file_size, rec.classification_reason))
                    if phase_group:
                        phase_rec = phase_group[echo_idx]
                        dst_phase_nii, _ = copy_phase_with_sidecar(phase_rec, dst_nii, args.dry_run, reporter)
                        manifest_rows.append(("sbref_phase", phase_rec.series_number, phase_rec.series_description, str(phase_rec.nii_path), str(phase_rec.json_path), str(dst_phase_nii), phase_rec.n_volumes, phase_rec.file_size, phase_rec.classification_reason))

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

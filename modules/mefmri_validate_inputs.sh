#!/usr/bin/env bash
# Input preflight validation for RevisedMe-fMRIPipeline.
# Call signature:
#   mefmri_validate_inputs.sh <SubjectDir> [FuncDirName] [FuncFilePrefix] [StartSession]

set -euo pipefail
IFS=$'\n\t'

SubjectDir="${1:?missing SubjectDir}"
FuncDirName="${2:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${3:-${FUNC_FILE_PREFIX:-Rest}}"
StartSession="${4:-${START_SESSION:-1}}"
VALIDATE_ECHO_DIM4_POLICY="${VALIDATE_ECHO_DIM4_POLICY:-error}"
FUNC_NOFIELDMAP_MODE="${FUNC_NOFIELDMAP_MODE:-0}"
DISTORTION_CORRECTION_MODE="${DISTORTION_CORRECTION_MODE:-topup}"
MEPCA="${MEPCA:-kundu}"

# A numeric tedana PCA setting is a requested component count. Only enforce
# that setting when tedana/ME-ICA is the effective denoising method. When this
# module is run directly, assume that a supplied MEPCA value is intended for
# tedana.
TEDANA_PCA_CHECK_ENABLED=1
if [[ -n "${PIPELINE_EFFECTIVE_DENOISE_METHOD:-}" && "${PIPELINE_EFFECTIVE_DENOISE_METHOD}" != "meica" ]]; then
  TEDANA_PCA_CHECK_ENABLED=0
fi
TEDANA_PCA_NCOMPS=""
if [[ "$TEDANA_PCA_CHECK_ENABLED" == "1" && "$MEPCA" =~ ^[0-9]+$ ]]; then
  TEDANA_PCA_NCOMPS="$((10#$MEPCA))"
fi

[[ -d "$SubjectDir" ]] || { echo "ERROR: missing subject directory: $SubjectDir" >&2; exit 2; }
[[ "$StartSession" =~ ^[0-9]+$ ]] || { echo "ERROR: StartSession must be integer, got: $StartSession" >&2; exit 2; }
[[ "$VALIDATE_ECHO_DIM4_POLICY" == "error" || "$VALIDATE_ECHO_DIM4_POLICY" == "warn" ]] || {
  echo "ERROR: VALIDATE_ECHO_DIM4_POLICY must be error|warn, got: $VALIDATE_ECHO_DIM4_POLICY" >&2
  exit 2
}
[[ "$FUNC_NOFIELDMAP_MODE" == "0" || "$FUNC_NOFIELDMAP_MODE" == "1" ]] || {
  echo "ERROR: FUNC_NOFIELDMAP_MODE must be 0|1, got: $FUNC_NOFIELDMAP_MODE" >&2
  exit 2
}

Subject="$(basename "$SubjectDir")"

ANAT_UNPROC="$SubjectDir/anat/unprocessed"
FUNC_UNPROC="$SubjectDir/func/unprocessed/$FuncDirName"
FM_UNPROC="${FM_RAW_DIR_REL:-func/unprocessed/${FuncDirName}/field_maps}"
FM_UNPROC="$SubjectDir/$FM_UNPROC"

declare -a ERRORS=()
declare -a WARNS=()

err() { ERRORS+=("$*"); }
warn() { WARNS+=("$*"); }

check_anat() {
  local t1_dir="$ANAT_UNPROC/T1w"
  local t2_dir="$ANAT_UNPROC/T2w"
  local t1_count t2_count
  [[ -d "$ANAT_UNPROC" ]] || { err "Missing anat raw directory: $ANAT_UNPROC"; return; }
  [[ -d "$t1_dir" ]] || { err "Missing T1w directory: $t1_dir"; return; }

  t1_count=$(find -L "$t1_dir" -maxdepth 1 -type f -name 'T1w*.nii.gz' | wc -l | tr -d ' ')
  if [[ "$t1_count" -lt 1 ]]; then
    err "No T1w*.nii.gz files found in $t1_dir"
  fi

  if [[ -d "$t2_dir" ]]; then
    t2_count=$(find -L "$t2_dir" -maxdepth 1 -type f -name 'T2w*.nii.gz' | wc -l | tr -d ' ')
    if [[ "$t2_count" -lt 1 ]]; then
      warn "T2w directory exists but no T2w*.nii.gz found: $t2_dir (pipeline can run in legacy mode)"
    fi
  else
    warn "Missing optional T2w directory: $t2_dir (pipeline will use legacy anatomical mode)"
  fi
}

check_func() {
  local session_dirs=() session run_dirs run run_num expected_prefix
  [[ -d "$FUNC_UNPROC" ]] || { err "Missing functional raw directory: $FUNC_UNPROC"; return; }

  mapfile -t session_dirs < <(find -L "$FUNC_UNPROC" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V)
  [[ "${#session_dirs[@]}" -gt 0 ]] || { err "No session_* folders found in $FUNC_UNPROC"; return; }

  for session in "${session_dirs[@]}"; do
    [[ "$session" =~ /session_([0-9]+)$ ]] || { err "Invalid session folder name (expected session_<N>): $session"; continue; }
    local s_num="${BASH_REMATCH[1]}"
    if (( s_num < StartSession )); then
      continue
    fi

    mapfile -t run_dirs < <(find -L "$session" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
    [[ "${#run_dirs[@]}" -gt 0 ]] || { err "No run_* folders found in $session"; continue; }

    for run in "${run_dirs[@]}"; do
      [[ "$run" =~ /run_([0-9]+)$ ]] || { err "Invalid run folder name (expected run_<N>): $run"; continue; }
      run_num="${BASH_REMATCH[1]}"
      expected_prefix="${FuncFilePrefix}_S${s_num}_R${run_num}_E"
      validate_run_payload "$run" "$expected_prefix"
    done
  done
}

validate_run_payload() {
  local run_dir="$1"
  local expected_prefix="$2"
  local -a nii_files=() json_files=()
  local nf jf echo_tag nvol
  local -a echo_sizes=() echo_names=() echo_vols=()
  local fslnvols_available=0

  if command -v fslnvols >/dev/null 2>&1; then
    fslnvols_available=1
  fi

  mapfile -t nii_files < <(find -L "$run_dir" -maxdepth 1 -type f -name "${FuncFilePrefix}_S*_R*_E*.nii.gz" | sed -E '/_phase\.nii\.gz$/d' | sort -V)
  if [[ "${#nii_files[@]}" -eq 0 ]]; then
    err "No echo NIfTI files found in $run_dir (expected ${expected_prefix}*.nii.gz)"
    return
  fi

  mapfile -t json_files < <(find -L "$run_dir" -maxdepth 1 -type f -name "${FuncFilePrefix}_S*_R*_E*.json" | sed -E '/_phase\.json$/d' | sort -V)
  if [[ "${#json_files[@]}" -eq 0 ]]; then
    err "No echo JSON sidecars found in $run_dir (expected ${expected_prefix}*.json)"
  fi

  for nf in "${nii_files[@]}"; do
    jf="${nf%.nii.gz}.json"
    if [[ ! -f "$jf" ]]; then
      err "Missing JSON sidecar for NIfTI: $nf"
    fi
    if [[ "$(basename "$nf")" != ${expected_prefix}* ]]; then
      err "Unexpected file naming in $run_dir: $(basename "$nf") (expected prefix ${expected_prefix})"
    fi
    echo_tag="$(basename "$nf" | sed -E 's/.*_E([0-9]+)\.nii\.gz/\1/')"
    if [[ ! "$echo_tag" =~ ^[0-9]+$ ]]; then
      err "Could not parse echo index from filename: $nf"
    fi
    echo_sizes+=("$(stat -Lc%s "$nf")")
    echo_names+=("$(basename "$nf")")
    if [[ "$fslnvols_available" == "1" ]]; then
      nvol=""
      if ! nvol="$(fslnvols "$nf" 2>/dev/null)"; then
        nvol=""
      fi
      echo_vols+=("$nvol")
    else
      echo_vols+=("")
    fi
  done

  if [[ "${#echo_sizes[@]}" -gt 1 ]]; then
    local last_idx=$(( ${#echo_sizes[@]} - 1 ))
    if (( echo_sizes[0] <= echo_sizes[last_idx] )); then
      warn "Echo size trend check in $run_dir: first echo (${echo_names[0]}, ${echo_sizes[0]} B) is not larger than last echo (${echo_names[last_idx]}, ${echo_sizes[last_idx]} B). Verify echo ordering/content."
    fi
    local i
    for (( i=0; i<last_idx; i++ )); do
      if (( echo_sizes[i] <= echo_sizes[i+1] )); then
        warn "Echo size trend check in $run_dir: expected decreasing size but found ${echo_names[i]} (${echo_sizes[i]} B) <= ${echo_names[i+1]} (${echo_sizes[i+1]} B). Double-check data ordering."
      fi
    done
  fi

  if [[ "${#echo_vols[@]}" -gt 0 ]]; then
    if [[ "$fslnvols_available" != "1" ]]; then
      warn "Echo dim4 and tedana PCA/TR checks skipped in $run_dir because fslnvols is unavailable on PATH."
    else
      local min_vol=""
      local max_vol=""
      local all_vols_numeric=1
      local required_trs=""
      local max_supported_components=""
      local i v

      if [[ -n "$TEDANA_PCA_NCOMPS" && "$TEDANA_PCA_NCOMPS" -gt 0 ]]; then
        # ceil(1.5 * n_components), calculated with integer arithmetic.
        required_trs="$(( (3 * TEDANA_PCA_NCOMPS + 1) / 2 ))"
      fi

      # Check every echo independently. A truncated final echo is common, but
      # any echo can be the offending series.
      for i in "${!echo_vols[@]}"; do
        v="${echo_vols[i]}"
        if [[ -z "$v" || ! "$v" =~ ^[0-9]+$ ]]; then
          warn "Could not determine dim4/TR count for ${echo_names[i]} in $run_dir (fslnvols output: '$v')."
          all_vols_numeric=0
          continue
        fi

        if [[ -z "$min_vol" || "$v" -lt "$min_vol" ]]; then
          min_vol="$v"
        fi
        if [[ -z "$max_vol" || "$v" -gt "$max_vol" ]]; then
          max_vol="$v"
        fi

        if [[ -n "$required_trs" && "$v" -lt "$required_trs" ]]; then
          max_supported_components="$(( (2 * v) / 3 ))"
          err "Tedana PCA/TR requirement failed for ${echo_names[i]} in $run_dir: found $v TRs, but numeric MEPCA=$TEDANA_PCA_NCOMPS requires at least $required_trs TRs (TRs must be at least 150% of the requested components). Use MEPCA<=$max_supported_components for this echo, choose a named tedana PCA method, or provide a longer functional series."
        fi
      done

      if [[ "$all_vols_numeric" == "1" && "${#echo_vols[@]}" -gt 1 && "$max_vol" -ne "$min_vol" ]]; then
        local details=""
        for i in "${!echo_vols[@]}"; do
          details+="${echo_names[i]}=${echo_vols[i]} "
        done
        local msg="Echo dim4 mismatch in $run_dir (per-echo timepoints: $details; min=$min_vol, max=$max_vol). This commonly indicates incomplete copy/truncation from source transfer. Recommended: re-import after fixing source data. Alternative: trim all echoes in this run to min=$min_vol (for example with fslroi) before running the pipeline."
        if [[ "$VALIDATE_ECHO_DIM4_POLICY" == "warn" ]]; then
          warn "$msg"
        else
          err "$msg"
        fi
      fi
    fi
  fi

  "${VALIDATE_PYTHON:-${PIPELINE_PYTHON:-python3}}" - "$run_dir" "$FuncFilePrefix" "$FUNC_NOFIELDMAP_MODE" <<'PY'
import json
import re
import sys
from pathlib import Path

run_dir = Path(sys.argv[1])
prefix = sys.argv[2]
no_fieldmap_mode = sys.argv[3] == "1"
pattern = re.compile(rf"^{re.escape(prefix)}_S\d+_R\d+_E\d+\.json$")
required = [
    "EchoTime",
    "RepetitionTime",
    "SliceTiming",
]
if not no_fieldmap_mode:
    required.extend(["EffectiveEchoSpacing", "PhaseEncodingDirection"])
errors = []

def parse_echo_idx(name: str):
    m = re.search(r"_E(\d+)\.(nii\.gz|json)$", name)
    return int(m.group(1)) if m else None

json_files = sorted([p for p in run_dir.glob(f"{prefix}_S*_R*_E*.json") if pattern.match(p.name)])
nii_pattern = re.compile(rf"^{re.escape(prefix)}_S\d+_R\d+_E\d+\.nii\.gz$")
nii_files = sorted([p for p in run_dir.glob(f"{prefix}_S*_R*_E*.nii.gz") if nii_pattern.match(p.name)])

nii_echoes = sorted([e for e in (parse_echo_idx(p.name) for p in nii_files) if e is not None])
json_echoes = sorted([e for e in (parse_echo_idx(p.name) for p in json_files) if e is not None])

if nii_echoes:
    expected = list(range(1, len(nii_echoes) + 1))
    if nii_echoes != expected:
        errors.append(
            f"Non-continuous echo numbering in {run_dir}: found {nii_echoes}, expected {expected}"
        )
if json_echoes:
    expected = list(range(1, len(json_echoes) + 1))
    if json_echoes != expected:
        errors.append(
            f"Non-continuous JSON echo numbering in {run_dir}: found {json_echoes}, expected {expected}"
        )
if nii_echoes and json_echoes and nii_echoes != json_echoes:
    errors.append(
        f"Echo mismatch between NIfTI and JSON in {run_dir}: nii={nii_echoes}, json={json_echoes}"
    )

for jf in json_files:
    try:
        data = json.loads(jf.read_text())
    except Exception as exc:
        errors.append(f"Invalid JSON ({jf}): {exc}")
        continue
    for key in required:
        if key not in data:
            errors.append(f"Missing key '{key}' in {jf}")
    if "SliceTiming" in data and not isinstance(data["SliceTiming"], list):
        errors.append(f"SliceTiming must be an array in {jf}")
    if "PhaseEncodingDirection" in data:
        ped = str(data["PhaseEncodingDirection"])
        if ped not in {"i", "i-", "j", "j-", "k", "k-"}:
            errors.append(f"Invalid PhaseEncodingDirection='{ped}' in {jf}")

if errors:
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(2)
PY
}

check_fieldmaps() {
  local -a ap_files=() pa_files=()
  local ap pa expected_json

  if [[ "$FUNC_NOFIELDMAP_MODE" == "1" || "$DISTORTION_CORRECTION_MODE" == "none" ]]; then
    if [[ -d "$FM_UNPROC" ]]; then
      warn "FUNC_NOFIELDMAP_MODE=1: fieldmap directory will be ignored: $FM_UNPROC"
    else
      warn "FUNC_NOFIELDMAP_MODE=1: skipping fieldmap presence checks."
    fi
    return
  fi

  if [[ "$DISTORTION_CORRECTION_MODE" == "medic" ]]; then
    warn "DISTORTION_CORRECTION_MODE=medic: skipping AP/PA fieldmap checks; MEDIC phase-companion checks are handled by the pipeline runner."
    return
  fi

  [[ -d "$FM_UNPROC" ]] || { err "Missing fieldmap raw directory: $FM_UNPROC"; return; }

  mapfile -t ap_files < <(find -L "$FM_UNPROC" -maxdepth 1 -type f -name 'AP_S*_R*.nii.gz' | sort -V)
  mapfile -t pa_files < <(find -L "$FM_UNPROC" -maxdepth 1 -type f -name 'PA_S*_R*.nii.gz' | sort -V)

  [[ "${#ap_files[@]}" -gt 0 ]] || err "No AP fieldmap NIfTI files found in $FM_UNPROC"
  [[ "${#pa_files[@]}" -gt 0 ]] || err "No PA fieldmap NIfTI files found in $FM_UNPROC"

  for ap in "${ap_files[@]}"; do
    pa="${ap/AP_/PA_}"
    [[ -f "$pa" ]] || err "Missing paired PA fieldmap for AP file: $ap"
    expected_json="${ap%.nii.gz}.json"
    [[ -f "$expected_json" ]] || err "Missing AP fieldmap JSON: $expected_json"
    expected_json="${pa%.nii.gz}.json"
    [[ -f "$expected_json" ]] || err "Missing PA fieldmap JSON: $expected_json"
  done
}

echo
echo "[validate] Subject: $Subject"
echo "[validate] SubjectDir: $SubjectDir"
echo "[validate] Functional naming: func/$FuncDirName, prefix ${FuncFilePrefix}_*"
echo "[validate] StartSession: $StartSession"
echo "[validate] Echo dim4 policy: $VALIDATE_ECHO_DIM4_POLICY"
echo "[validate] Distortion correction mode: $DISTORTION_CORRECTION_MODE"
echo "[validate] FUNC_NOFIELDMAP_MODE: $FUNC_NOFIELDMAP_MODE"
if [[ "$TEDANA_PCA_CHECK_ENABLED" == "1" ]]; then
  echo "[validate] Tedana PCA method: $MEPCA"
else
  echo "[validate] Tedana PCA/TR check: not applicable (effective denoising method: ${PIPELINE_EFFECTIVE_DENOISE_METHOD})"
fi

if [[ -n "$TEDANA_PCA_NCOMPS" && "$TEDANA_PCA_NCOMPS" -lt 1 ]]; then
  err "Numeric MEPCA must request at least one component, got: $MEPCA"
fi

check_anat
check_func
check_fieldmaps

if [[ "${#WARNS[@]}" -gt 0 ]]; then
  echo
  echo "[validate] WARNINGS"
  for w in "${WARNS[@]}"; do
    echo "WARNING: $w"
  done
fi

if [[ "${#ERRORS[@]}" -gt 0 ]]; then
  echo
  echo "[validate] ERRORS (${#ERRORS[@]})"
  for e in "${ERRORS[@]}"; do
    echo "ERROR: $e"
  done
  echo
  echo "[validate] Input validation failed. Fix the errors above and re-run."
  exit 2
fi

echo
echo "[validate] Input validation passed."

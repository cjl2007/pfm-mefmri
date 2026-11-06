#!/usr/bin/env bash
# Unified importer front-end for raw DICOM exports and BIDS NIfTI datasets.
#
# Raw DICOM inputs are dispatched to mefmri_import_raw.sh.
# BIDS/NIfTI inputs are dispatched to mefmri_import_bids.sh.

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  mefmri_import.sh [--input-type auto|raw|bids] <INPUT> <...backend args...>

Common examples:
  # Raw scanner export/DICOM directory. Existing raw-import arguments are passed through.
  mefmri_import.sh /path/to/raw_dicom_export /path/to/study/ME001 config/mefmri_import_raw_config.sh --session 1 --nordic

  # Multiple raw visits. Existing raw-import arguments are passed through.
  mefmri_import.sh --raw-parent /path/to/Raw --subject-dir /path/to/study/ME001 config/mefmri_import_raw_config.sh --nordic

  # BIDS root using legacy BIDS-import positional arguments.
  mefmri_import.sh /path/to/bids 06 /path/to/study/ME06 --task rest --mode symlink --overwrite

  # BIDS root or sub-<label> directory using named arguments.
  mefmri_import.sh /path/to/bids --subject 06 --subject-dir /path/to/study/ME06 --task rest --mode symlink --overwrite
  mefmri_import.sh /path/to/bids/sub-06 --subject-dir /path/to/study/ME06 --task rest --mode symlink --overwrite

Input detection:
  auto: raw if DICOM-like files (*.dcm, DICOMDIR) are found; BIDS if BIDS-style
        sub-* NIfTIs are found.
  raw:  force raw importer.
  bids: force BIDS importer.

NORDIC:
  --nordic is forwarded to the raw importer and runs/stages NORDIC there.
  For BIDS imports, NORDIC metadata already present in BIDS sidecars is preserved
  and converted into per-run NORDIC_DENOISING.txt markers. Running NORDIC from
  BIDS NIfTIs is future work.
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_IMPORTER="$SCRIPT_DIR/mefmri_import_raw.sh"
BIDS_IMPORTER="$SCRIPT_DIR/mefmri_import_bids.sh"

die() {
  echo "ERROR: $*" >&2
  exit 2
}

strip_unified_args() {
  local -n _in="$1"
  local -n _out="$2"
  local i=0
  _out=()
  while (( i < ${#_in[@]} )); do
    case "${_in[$i]}" in
      --input-type|--type)
        (( i + 1 < ${#_in[@]} )) || die "${_in[$i]} requires auto|raw|bids"
        i=$((i + 2))
        ;;
      --input-type=*|--type=*)
        i=$((i + 1))
        ;;
      *)
        _out+=("${_in[$i]}")
        i=$((i + 1))
        ;;
    esac
  done
}

first_positional() {
  local args=("$@")
  local i=0
  while (( i < ${#args[@]} )); do
    case "${args[$i]}" in
      --input-type|--type|--subject|--bids-subject|--subject-dir|--raw-dir|--raw-list|--raw-parent|--session|--import-work-subdir|--nordic-code-dir|--nordic-matlab-bin|--task|--func-dirname|--func-prefix|--mode|--echo-dim4-policy|--expect-task-volumes)
        i=$((i + 2))
        ;;
      --*=*|--dry-run|--nordic|-nordic|--prepare-nordic-commands-only|--make-nordic-test-sessions|--overwrite)
        i=$((i + 1))
        ;;
      --*)
        i=$((i + 1))
        ;;
      *)
        printf '%s\n' "${args[$i]}"
        return 0
        ;;
    esac
  done
  return 1
}

contains_dicoms() {
  local p="$1"
  [[ -e "$p" ]] || return 1
  if [[ -f "$p" ]]; then
    [[ "${p,,}" == *.dcm || "$(basename "$p")" == "DICOMDIR" ]]
    return
  fi
  find "$p" -maxdepth 4 \( -iname '*.dcm' -o -name 'DICOMDIR' \) -print -quit 2>/dev/null | grep -q .
}

looks_like_bids() {
  local p="$1"
  [[ -e "$p" ]] || return 1
  if [[ -f "$p" ]]; then
    [[ "${p,,}" == *.nii || "${p,,}" == *.nii.gz ]] || return 1
    return 1
  fi
  if [[ -f "$p/dataset_description.json" ]]; then
    return 0
  fi
  if [[ "$(basename "$p")" == sub-* && ( -d "$p/func" || -d "$p/anat" ) ]]; then
    find "$p" -maxdepth 3 \( -iname '*.nii' -o -iname '*.nii.gz' \) -print -quit 2>/dev/null | grep -q .
    return
  fi
  find "$p" -maxdepth 4 -path '*/sub-*/*' \( -iname '*.nii' -o -iname '*.nii.gz' \) -print -quit 2>/dev/null | grep -q .
}

get_arg_value() {
  local key="$1"
  shift
  local args=("$@")
  local i=0
  while (( i < ${#args[@]} )); do
    case "${args[$i]}" in
      "$key")
        (( i + 1 < ${#args[@]} )) || die "$key requires a value"
        printf '%s\n' "${args[$((i + 1))]}"
        return 0
        ;;
      "$key="*)
        printf '%s\n' "${args[$i]#*=}"
        return 0
        ;;
    esac
    i=$((i + 1))
  done
  return 1
}

has_flag() {
  local key="$1"
  shift
  local a
  for a in "$@"; do
    [[ "$a" == "$key" ]] && return 0
  done
  return 1
}

build_bids_args() {
  local -n _in="$1"
  local -n _out="$2"
  local pos=()
  local rest=()
  local subject=""
  local subject_dir=""
  local saw_nordic=0
  local i=0

  _out=()
  while (( i < ${#_in[@]} )); do
    case "${_in[$i]}" in
      --input-type|--type)
        i=$((i + 2))
        ;;
      --input-type=*|--type=*)
        i=$((i + 1))
        ;;
      --subject|--bids-subject)
        (( i + 1 < ${#_in[@]} )) || die "${_in[$i]} requires a value"
        subject="${_in[$((i + 1))]}"
        i=$((i + 2))
        ;;
      --subject=*|--bids-subject=*)
        subject="${_in[$i]#*=}"
        i=$((i + 1))
        ;;
      --subject-dir)
        (( i + 1 < ${#_in[@]} )) || die "--subject-dir requires a value"
        subject_dir="${_in[$((i + 1))]}"
        i=$((i + 2))
        ;;
      --subject-dir=*)
        subject_dir="${_in[$i]#*=}"
        i=$((i + 1))
        ;;
      --nordic|-nordic)
        saw_nordic=1
        i=$((i + 1))
        ;;
      --prepare-nordic-commands-only|--make-nordic-test-sessions|--nordic-code-dir|--nordic-matlab-bin)
        die "${_in[$i]} is raw-import only; BIDS import preserves existing NORDIC metadata but does not run NORDIC."
        ;;
      --*)
        rest+=("${_in[$i]}")
        if [[ "${_in[$i]}" != *=* && "$((i + 1))" -lt "${#_in[@]}" && "${_in[$((i + 1))]}" != --* ]]; then
          case "${_in[$i]}" in
            --task|--func-dirname|--func-prefix|--mode|--echo-dim4-policy|--expect-task-volumes)
              rest+=("${_in[$((i + 1))]}")
              i=$((i + 2))
              continue
              ;;
          esac
        fi
        i=$((i + 1))
        ;;
      *)
        pos+=("${_in[$i]}")
        i=$((i + 1))
        ;;
    esac
  done

  if (( saw_nordic )); then
    echo "[import] BIDS --nordic requested: preserving existing NORDIC metadata only; BIDS-side NORDIC execution is not implemented." >&2
  fi

  if (( ${#pos[@]} >= 3 )) && [[ -z "$subject" && -z "$subject_dir" ]]; then
    _out=("${pos[0]}" "${pos[1]}" "${pos[2]}" "${rest[@]}")
    return 0
  fi

  (( ${#pos[@]} >= 1 )) || die "BIDS import requires BIDS_ROOT or sub-<subject> directory"
  local bids_input="${pos[0]}"
  if [[ "$(basename "$bids_input")" == sub-* && ( -d "$bids_input/func" || -d "$bids_input/anat" ) ]]; then
    subject="${subject:-$(basename "$bids_input")}"
    bids_input="$(cd "$bids_input/.." && pwd)"
  fi
  [[ -n "$subject" ]] || die "BIDS import requires --subject LABEL, or use a sub-<label> input directory, or use legacy <BIDS_ROOT> <SUBJECT> <OUT_SUBJECT_DIR> syntax"
  [[ -n "$subject_dir" ]] || die "BIDS import requires --subject-dir OUT_SUBJECT_DIR, or use legacy <BIDS_ROOT> <SUBJECT> <OUT_SUBJECT_DIR> syntax"
  _out=("$bids_input" "$subject" "$subject_dir" "${rest[@]}")
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -lt 1 ]]; then
  usage >&2
  exit 2
fi

INPUT_TYPE="auto"
case "${1:-}" in
  --input-type|--type)
    [[ "$#" -ge 2 ]] || die "$1 requires auto|raw|bids"
    INPUT_TYPE="$2"
    ;;
  --input-type=*|--type=*)
    INPUT_TYPE="${1#*=}"
    ;;
esac
case "$INPUT_TYPE" in
  auto|raw|bids) ;;
  *) die "--input-type must be auto, raw, or bids (got '$INPUT_TYPE')" ;;
esac

ARGS=("$@")
DISPATCH_ARGS=()

if [[ "$INPUT_TYPE" == "auto" ]]; then
  if has_flag "--raw-dir" "${ARGS[@]}" || has_flag "--raw-list" "${ARGS[@]}" || has_flag "--raw-parent" "${ARGS[@]}"; then
    INPUT_TYPE="raw"
  else
    input="$(first_positional "${ARGS[@]}" || true)"
    [[ -n "${input:-}" ]] || die "Could not infer input type without an input path"
    if contains_dicoms "$input"; then
      INPUT_TYPE="raw"
    elif looks_like_bids "$input"; then
      INPUT_TYPE="bids"
    else
      die "Could not infer input type for '$input'. Use --input-type raw or --input-type bids."
    fi
  fi
fi

case "$INPUT_TYPE" in
  raw)
    strip_unified_args ARGS DISPATCH_ARGS
    echo "[import] detected input_type=raw"
    exec "$RAW_IMPORTER" "${DISPATCH_ARGS[@]}"
    ;;
  bids)
    build_bids_args ARGS DISPATCH_ARGS
    echo "[import] detected input_type=bids"
    exec "$BIDS_IMPORTER" "${DISPATCH_ARGS[@]}"
    ;;
esac

#!/usr/bin/env bash
# Convert raw scanner-export DICOM folders into RevisedMe-fMRIPipeline raw layout.
#
# Usage:
#   mefmri_import_raw.sh <RAW_DICOM_DIR> <SUBJECT_DIR> [CONFIG_FILE] [--session N] [--dry-run]
#
# Example:
#   bash bin/mefmri_import_raw.sh \
#     /path/to/raw_dicom_export \
#     /path/to/study/ME001 \
#     /path/to/config/mefmri_import_raw_bd2_config.sh \
#     --session 1
#
# Dry-run example:
#   bash bin/mefmri_import_raw.sh \
#     /path/to/raw_dicom_export \
#     /path/to/study/ME001 \
#     /path/to/config/mefmri_import_raw_bd2_config.sh \
#     --session 1 \
#     --dry-run

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  mefmri_import_raw.sh <RAW_DICOM_DIR> <SUBJECT_DIR> [CONFIG_FILE] [--session N] [--dry-run]

Example:
  bash bin/mefmri_import_raw.sh \
    /path/to/raw_dicom_export \
    /path/to/study/ME001 \
    config/mefmri_import_raw_bd2_config.sh \
    --dry-run
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -lt 2 ]]; then
  usage >&2
  exit 2
fi

RAW_DICOM_DIR="$1"
SUBJECT_DIR="$2"
shift 2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$MEDIR/config/mefmri_import_config_current_protocol.sh"

if [[ "${1:-}" != "" && "${1:-}" != --* ]]; then
  CONFIG_FILE="$1"
  shift
fi

SESSION=""
DRY_RUN=0

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --session)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --session requires an integer value." >&2; exit 2; }
      SESSION="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -d "$RAW_DICOM_DIR" ]] || { echo "ERROR: missing raw DICOM directory: $RAW_DICOM_DIR" >&2; exit 2; }
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: missing config file: $CONFIG_FILE" >&2; exit 2; }

source "$CONFIG_FILE"

serialize_array() {
  local array_name="$1"
  local value=""
  if declare -p "$array_name" >/dev/null 2>&1; then
    local -n arr_ref="$array_name"
    local item
    for item in "${arr_ref[@]}"; do
      if [[ -n "$value" ]]; then
        value+=$'\n'
      fi
      value+="$item"
    done
  fi
  printf '%s' "$value"
}

next_session_suggestion() {
  local subject_dir="$1"
  local max_session=0
  local session_dirs=()
  if [[ -d "$subject_dir/func/unprocessed" ]]; then
    mapfile -t session_dirs < <(find "$subject_dir/func/unprocessed" -mindepth 2 -maxdepth 2 -type d -name 'session_*' | sort -V)
    local sdir sval
    for sdir in "${session_dirs[@]}"; do
      sval="${sdir##*/}"
      sval="${sval#session_}"
      if [[ "$sval" =~ ^[0-9]+$ ]] && (( sval > max_session )); then
        max_session="$sval"
      fi
    done
  fi
  printf '%s' "$((max_session + 1))"
}

if [[ -n "$SESSION" ]] && [[ ! "$SESSION" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --session must be an integer, got: $SESSION" >&2
  exit 2
fi

if [[ -d "$SUBJECT_DIR" && -z "$SESSION" ]]; then
  suggested_session="$(next_session_suggestion "$SUBJECT_DIR")"
  if [[ -t 0 ]]; then
    read -r -p "Subject dir exists. Enter session number [${suggested_session}]: " reply
    SESSION="${reply:-$suggested_session}"
  else
    echo "ERROR: subject dir exists: $SUBJECT_DIR" >&2
    echo "Re-run with --session ${suggested_session}" >&2
    exit 2
  fi
fi

if [[ -z "$SESSION" ]]; then
  SESSION=1
fi

[[ "$SESSION" =~ ^[0-9]+$ ]] || { echo "ERROR: session must be an integer, got: $SESSION" >&2; exit 2; }

export IMPORT_PROTOCOL_NAME="${IMPORT_PROTOCOL_NAME:-}"
export IMPORT_DCM2NIIX_BIN="${IMPORT_DCM2NIIX_BIN:-dcm2niix}"
export FUNC_DIRNAME="${FUNC_DIRNAME:-rest}"
export FUNC_FILE_PREFIX="${FUNC_FILE_PREFIX:-Rest}"
export IMPORT_EXPECT_REST_RUNS_PER_SESSION="${IMPORT_EXPECT_REST_RUNS_PER_SESSION:-0}"
export IMPORT_EXPECT_ECHOES_PER_RUN="${IMPORT_EXPECT_ECHOES_PER_RUN:-0}"
export IMPORT_EXPECT_SBREF_PER_RUN="${IMPORT_EXPECT_SBREF_PER_RUN:-0}"
export IMPORT_EXPECT_FMAP_AP_PER_SESSION="${IMPORT_EXPECT_FMAP_AP_PER_SESSION:-0}"
export IMPORT_EXPECT_FMAP_PA_PER_SESSION="${IMPORT_EXPECT_FMAP_PA_PER_SESSION:-0}"
export IMPORT_EXPECT_T1W_MAX_PER_IMPORT="${IMPORT_EXPECT_T1W_MAX_PER_IMPORT:-0}"
export IMPORT_EXPECT_T2W_MAX_PER_IMPORT="${IMPORT_EXPECT_T2W_MAX_PER_IMPORT:-0}"
export IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING="${IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING:-0}"
export IMPORT_T2W_OPTIONAL="${IMPORT_T2W_OPTIONAL:-1}"
export IMPORT_EXPECT_REST_VOLUMES="${IMPORT_EXPECT_REST_VOLUMES:-0}"
export IMPORT_EXPECT_SBREF_VOLUMES="${IMPORT_EXPECT_SBREF_VOLUMES:-0}"
export IMPORT_EXPECT_FMAP_VOLUMES="${IMPORT_EXPECT_FMAP_VOLUMES:-0}"
export IMPORT_EXPECT_ANAT_VOLUMES="${IMPORT_EXPECT_ANAT_VOLUMES:-0}"
export IMPORT_MIN_BYTES_REST="${IMPORT_MIN_BYTES_REST:-0}"
export IMPORT_MIN_BYTES_SBREF="${IMPORT_MIN_BYTES_SBREF:-0}"
export IMPORT_MIN_BYTES_FMAP="${IMPORT_MIN_BYTES_FMAP:-0}"
export IMPORT_MIN_BYTES_T1W="${IMPORT_MIN_BYTES_T1W:-0}"
export IMPORT_MIN_BYTES_T2W="${IMPORT_MIN_BYTES_T2W:-0}"
export IMPORT_T1W_REGEX="${IMPORT_T1W_REGEX:-}"
export IMPORT_T2W_REGEX="${IMPORT_T2W_REGEX:-}"
export IMPORT_REST_REGEX="${IMPORT_REST_REGEX:-}"
export IMPORT_SBREF_REGEX="${IMPORT_SBREF_REGEX:-}"
export IMPORT_FMAP_AP_REGEX="${IMPORT_FMAP_AP_REGEX:-}"
export IMPORT_FMAP_PA_REGEX="${IMPORT_FMAP_PA_REGEX:-}"
export IMPORT_IGNORE_REGEXES_SERIALIZED="$(serialize_array IMPORT_IGNORE_REGEXES)"

python3 "$MEDIR/lib/mefmri_raw_import.py" \
  "$RAW_DICOM_DIR" \
  "$SUBJECT_DIR" \
  --session "$SESSION" \
  --config-file "$CONFIG_FILE" \
  $( [[ "$DRY_RUN" -eq 1 ]] && printf '%s' '--dry-run' )

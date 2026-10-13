#!/usr/bin/env bash
# Convert raw scanner-export DICOM folders into RevisedMe-fMRIPipeline raw layout.
#
# Usage:
#   mefmri_import_raw.sh <RAW_DICOM_DIR> <SUBJECT_DIR> [CONFIG_FILE] [--session N] [--import-work-subdir DIR] [--dry-run] [--nordic] [--prepare-nordic-commands-only]
#   mefmri_import_raw.sh --raw-dir DIR [--raw-dir DIR ...] --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]
#   mefmri_import_raw.sh --raw-list FILE --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]
#   mefmri_import_raw.sh --raw-parent DIR --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]
#
# Example (copy and edit config/mefmri_import_config.sh first):
#   bash bin/mefmri_import_raw.sh \
#     /path/to/raw_dicom_export \
#     /path/to/study/ME001 \
#     /path/to/config/mefmri_import_config.sh \
#     --session 1
#
# Dry-run example:
#   bash bin/mefmri_import_raw.sh \
#     /path/to/raw_dicom_export \
#     /path/to/study/ME001 \
#     /path/to/config/mefmri_import_config.sh \
#     --session 1 \
#     --dry-run

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  mefmri_import_raw.sh <RAW_DICOM_DIR> <SUBJECT_DIR> [CONFIG_FILE] [--session N] [--import-work-subdir DIR] [--dry-run] [--nordic] [--prepare-nordic-commands-only]
  mefmri_import_raw.sh --raw-dir DIR [--raw-dir DIR ...] --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]
  mefmri_import_raw.sh --raw-list FILE --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]
  mefmri_import_raw.sh --raw-parent DIR --subject-dir SUBJECT_DIR [CONFIG_FILE] [--session N] [--dry-run] [--nordic]

Example (copy and edit config/mefmri_import_config.sh first):
  bash bin/mefmri_import_raw.sh \
    /path/to/raw_dicom_export \
    /path/to/study/ME001 \
    /path/to/study/config/import.sh \
    --dry-run

Multiple visits:
  bash bin/mefmri_import_raw.sh \
    --raw-parent "/path/to/raw/SUB001" \
    --subject-dir /path/to/study/SUB001 \
    /path/to/study/config/import.sh \
    --nordic
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -lt 1 ]]; then
  usage >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$MEDIR/config/mefmri_import_config.sh"

SESSION=""
DRY_RUN=0
USE_NORDIC=0
MAKE_NORDIC_TEST_SESSIONS=0
PREPARE_NORDIC_COMMANDS_ONLY=0
IMPORT_WORK_SUBDIR=""
NORDIC_CODE_DIR="${NORDIC_CODE_DIR:-$MEDIR/lib/nordic}"
NORDIC_MATLAB_BIN="${NORDIC_MATLAB_BIN:-matlab}"
SUBJECT_DIR=""
RAW_PARENT=""
RAW_LIST_FILE=""
RAW_DICOM_DIRS=()

extract_date_sort_key() {
  local path="$1"
  local name
  name="$(basename "$path")"
  if [[ "$name" =~ ([12][0-9]{3})[-_]?([01][0-9])[-_]?([0-3][0-9]) ]]; then
    printf '%s%s%s\t%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" "${BASH_REMATCH[3]}" "$path"
  elif [[ "$name" =~ ([01][0-9])[-_]?([0-3][0-9])[-_]?([12][0-9]{3}) ]]; then
    printf '%s%s%s\t%s\n' "${BASH_REMATCH[3]}" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" "$path"
  else
    printf '99999999\t%s\n' "$path"
  fi
}

read_raw_list_file() {
  local list_file="$1"
  [[ -f "$list_file" ]] || { echo "ERROR: missing --raw-list file: $list_file" >&2; exit 2; }
  local line
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -n "$line" ]] || continue
    RAW_DICOM_DIRS+=("$line")
  done < "$list_file"
}

discover_raw_parent_dirs() {
  local parent="$1"
  [[ -d "$parent" ]] || { echo "ERROR: missing --raw-parent directory: $parent" >&2; exit 2; }
  local discovered=()
  mapfile -t discovered < <(
    find "$parent" -mindepth 1 -maxdepth 1 -type d -print0 |
      while IFS= read -r -d '' d; do
        extract_date_sort_key "$d"
      done |
      sort -k1,1 -k2,2 |
      cut -f2-
  )
  [[ "${#discovered[@]}" -gt 0 ]] || { echo "ERROR: no child directories found under --raw-parent: $parent" >&2; exit 2; }
  RAW_DICOM_DIRS+=("${discovered[@]}")
}

if [[ "${1:-}" != --* && "${1:-}" != "-nordic" ]]; then
  [[ "$#" -ge 2 ]] || { usage >&2; exit 2; }
  RAW_DICOM_DIRS+=("$1")
  SUBJECT_DIR="$2"
  shift 2
  if [[ "${1:-}" != "" && "${1:-}" != --* && "${1:-}" != "-nordic" ]]; then
    CONFIG_FILE="$1"
    shift
  fi
fi

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --raw-dir)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --raw-dir requires a path." >&2; exit 2; }
      RAW_DICOM_DIRS+=("$2")
      shift 2
      ;;
    --raw-list)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --raw-list requires a path." >&2; exit 2; }
      RAW_LIST_FILE="$2"
      shift 2
      ;;
    --raw-parent)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --raw-parent requires a path." >&2; exit 2; }
      RAW_PARENT="$2"
      shift 2
      ;;
    --subject-dir)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --subject-dir requires a path." >&2; exit 2; }
      SUBJECT_DIR="$2"
      shift 2
      ;;
    --session)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --session requires an integer value." >&2; exit 2; }
      SESSION="$2"
      shift 2
      ;;
    --import-work-subdir)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --import-work-subdir requires a relative path." >&2; exit 2; }
      IMPORT_WORK_SUBDIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --nordic|-nordic)
      USE_NORDIC=1
      shift
      ;;
    --make-nordic-test-sessions)
      echo "WARNING: --make-nordic-test-sessions is deprecated; use --nordic for normal imports. This debug flag still creates session_1=NORDIC and session_2=original for one raw directory." >&2
      MAKE_NORDIC_TEST_SESSIONS=1
      USE_NORDIC=1
      shift
      ;;
    --prepare-nordic-commands-only)
      PREPARE_NORDIC_COMMANDS_ONLY=1
      USE_NORDIC=1
      shift
      ;;
    --nordic-code-dir)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --nordic-code-dir requires a path." >&2; exit 2; }
      NORDIC_CODE_DIR="$2"
      shift 2
      ;;
    --nordic-matlab-bin)
      [[ "$#" -ge 2 ]] || { echo "ERROR: --nordic-matlab-bin requires a command/path." >&2; exit 2; }
      NORDIC_MATLAB_BIN="$2"
      shift 2
      ;;
    --*)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      if [[ -f "$1" && "$CONFIG_FILE" == "$MEDIR/config/mefmri_import_config.sh" ]]; then
        CONFIG_FILE="$1"
        shift
      else
        echo "ERROR: unexpected positional argument: $1" >&2
        usage >&2
        exit 2
      fi
      ;;
  esac
done

if [[ -n "$RAW_LIST_FILE" ]]; then
  read_raw_list_file "$RAW_LIST_FILE"
fi
if [[ -n "$RAW_PARENT" ]]; then
  discover_raw_parent_dirs "$RAW_PARENT"
fi

[[ "${#RAW_DICOM_DIRS[@]}" -gt 0 ]] || { echo "ERROR: no raw DICOM directories were provided." >&2; usage >&2; exit 2; }
[[ -n "$SUBJECT_DIR" ]] || { echo "ERROR: missing subject directory. Use positional SUBJECT_DIR or --subject-dir." >&2; usage >&2; exit 2; }
for raw_dir in "${RAW_DICOM_DIRS[@]}"; do
  [[ -d "$raw_dir" ]] || { echo "ERROR: missing raw DICOM directory: $raw_dir" >&2; exit 2; }
done
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

if [[ "$MAKE_NORDIC_TEST_SESSIONS" -eq 1 ]]; then
  [[ "${#RAW_DICOM_DIRS[@]}" -eq 1 ]] || { echo "ERROR: --make-nordic-test-sessions only supports one raw DICOM directory. Use --nordic for general multi-session imports." >&2; exit 2; }
  SESSION=1
fi

export IMPORT_PROTOCOL_NAME="${IMPORT_PROTOCOL_NAME:-}"
export IMPORT_DCM2NIIX_BIN="${IMPORT_DCM2NIIX_BIN:-dcm2niix}"
export FUNC_DIRNAME="${FUNC_DIRNAME:-rest}"
export FUNC_FILE_PREFIX="${FUNC_FILE_PREFIX:-Rest}"
export IMPORT_EXPECT_TASK_RUNS_PER_SESSION="${IMPORT_EXPECT_TASK_RUNS_PER_SESSION:-${IMPORT_EXPECT_REST_RUNS_PER_SESSION:-0}}"
export IMPORT_EXPECT_REST_RUNS_PER_SESSION="${IMPORT_EXPECT_TASK_RUNS_PER_SESSION}"
export IMPORT_EXPECT_ECHOES_PER_RUN="${IMPORT_EXPECT_ECHOES_PER_RUN:-0}"
export IMPORT_EXPECT_SBREF_PER_RUN="${IMPORT_EXPECT_SBREF_PER_RUN:-0}"
export IMPORT_EXPECT_FMAP_AP_PER_SESSION="${IMPORT_EXPECT_FMAP_AP_PER_SESSION:-0}"
export IMPORT_EXPECT_FMAP_PA_PER_SESSION="${IMPORT_EXPECT_FMAP_PA_PER_SESSION:-0}"
export IMPORT_EXPECT_T1W_MAX_PER_IMPORT="${IMPORT_EXPECT_T1W_MAX_PER_IMPORT:-0}"
export IMPORT_EXPECT_T2W_MAX_PER_IMPORT="${IMPORT_EXPECT_T2W_MAX_PER_IMPORT:-0}"
export IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING="${IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING:-0}"
export IMPORT_T2W_OPTIONAL="${IMPORT_T2W_OPTIONAL:-1}"
export IMPORT_EXPECT_TASK_VOLUMES="${IMPORT_EXPECT_TASK_VOLUMES:-${IMPORT_EXPECT_REST_VOLUMES:-0}}"
export IMPORT_EXPECT_REST_VOLUMES="${IMPORT_EXPECT_TASK_VOLUMES}"
export IMPORT_EXPECT_SBREF_VOLUMES="${IMPORT_EXPECT_SBREF_VOLUMES:-0}"
export IMPORT_EXPECT_FMAP_VOLUMES="${IMPORT_EXPECT_FMAP_VOLUMES:-0}"
export IMPORT_EXPECT_ANAT_VOLUMES="${IMPORT_EXPECT_ANAT_VOLUMES:-0}"
export IMPORT_ECHO_DIM4_POLICY="${IMPORT_ECHO_DIM4_POLICY:-abort}"
export IMPORT_MIN_BYTES_TASK="${IMPORT_MIN_BYTES_TASK:-${IMPORT_MIN_BYTES_REST:-0}}"
export IMPORT_MIN_BYTES_REST="${IMPORT_MIN_BYTES_TASK}"
export IMPORT_MIN_BYTES_SBREF="${IMPORT_MIN_BYTES_SBREF:-0}"
export IMPORT_MIN_BYTES_FMAP="${IMPORT_MIN_BYTES_FMAP:-0}"
export IMPORT_MIN_BYTES_T1W="${IMPORT_MIN_BYTES_T1W:-0}"
export IMPORT_MIN_BYTES_T2W="${IMPORT_MIN_BYTES_T2W:-0}"
export IMPORT_T1W_REGEX="${IMPORT_T1W_REGEX:-}"
export IMPORT_T2W_REGEX="${IMPORT_T2W_REGEX:-}"
export IMPORT_TASK_REGEX="${IMPORT_TASK_REGEX:-${IMPORT_REST_REGEX:-}}"
export IMPORT_REST_REGEX="${IMPORT_TASK_REGEX}"
export IMPORT_SBREF_REGEX="${IMPORT_SBREF_REGEX:-}"
export IMPORT_FMAP_AP_REGEX="${IMPORT_FMAP_AP_REGEX:-}"
export IMPORT_FMAP_PA_REGEX="${IMPORT_FMAP_PA_REGEX:-}"
export IMPORT_IGNORE_REGEXES_SERIALIZED="$(serialize_array IMPORT_IGNORE_REGEXES)"
PIPELINE_PYTHON="${PIPELINE_PYTHON:-python3}"

for idx in "${!RAW_DICOM_DIRS[@]}"; do
  raw_dir="${RAW_DICOM_DIRS[$idx]}"
  import_session="$((SESSION + idx))"
  IMPORT_ARGS=(
    "$raw_dir"
    "$SUBJECT_DIR"
    --session "$import_session"
    --config-file "$CONFIG_FILE"
  )
  if [[ "${#RAW_DICOM_DIRS[@]}" -gt 1 ]]; then
    IMPORT_ARGS+=(--import-work-subdir "session_${import_session}")
  elif [[ -n "$IMPORT_WORK_SUBDIR" ]]; then
    IMPORT_ARGS+=(--import-work-subdir "$IMPORT_WORK_SUBDIR")
  fi
  if [[ "$USE_NORDIC" -eq 1 ]]; then
    IMPORT_ARGS+=(--nordic)
    IMPORT_ARGS+=(--nordic-code-dir "$NORDIC_CODE_DIR")
    IMPORT_ARGS+=(--nordic-matlab-bin "$NORDIC_MATLAB_BIN")
    if [[ "$PREPARE_NORDIC_COMMANDS_ONLY" -eq 1 ]]; then
      IMPORT_ARGS+=(--prepare-nordic-commands-only)
    fi
  fi
  if [[ "$MAKE_NORDIC_TEST_SESSIONS" -eq 1 ]]; then
    IMPORT_ARGS+=(--make-nordic-test-sessions)
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    IMPORT_ARGS+=(--dry-run)
  fi

  echo "[import-wrapper] raw dir $((idx + 1))/${#RAW_DICOM_DIRS[@]} -> session_${import_session}: $raw_dir"
  "$PIPELINE_PYTHON" "$MEDIR/lib/mefmri_raw_import.py" "${IMPORT_ARGS[@]}"
done

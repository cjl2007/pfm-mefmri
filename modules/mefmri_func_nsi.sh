#!/usr/bin/env bash
# NSI module: run external NSI CLI from concatenated + censored CIFTI.
# Supports separate usability and reliability passes.
set -euo pipefail

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
MEDIR="${3:?missing MEDIR}"
StartSession="${4:?missing StartSession}"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

_unused_start_session="${StartSession}" # reserved for future non-concat NSI entrypoints
_unused_medir="${MEDIR}" # reserved for future in-tree NSI backend

NSI_ENABLE="${NSI_ENABLE:-${PFM_NSI_ENABLE:-1}}"
NSI_USE_EXTERNAL_CLI="${NSI_USE_EXTERNAL_CLI:-${PFM_USE_EXTERNAL_CLI:-1}}"
NSI_PYTHON="${NSI_PYTHON:-${PFM_PYTHON:-${PIPELINE_PYTHON:-python3}}}"
NSI_INPUT_TAG="${NSI_INPUT_TAG:-${CONCAT_INPUT_TAG:-${PFM_INPUT_TAG:-OCME+MEICA+MGTR}}}"
NSI_CONCAT_OUT_SUBDIR="${NSI_CONCAT_OUT_SUBDIR:-${CONCAT_OUT_SUBDIR:-${PFM_OUT_SUBDIR:-ConcatenatedCiftis}}}"
NSI_FD_THRESHOLD="${NSI_FD_THRESHOLD:-${CONCAT_FD_THRESHOLD:-${PFM_FD_THRESHOLD:-0.3}}}"

NSI_USABILITY_MODEL="${NSI_USABILITY_MODEL:-${NSI_EXTERNAL_USABILITY:-${PFM_EXTERNAL_USABILITY:-1}}}"
NSI_RELIABILITY_MODEL="${NSI_RELIABILITY_MODEL:-${NSI_EXTERNAL_RELIABILITY:-${PFM_EXTERNAL_RELIABILITY:-0}}}"
NSI_RELIABILITY_NSI_T="${NSI_RELIABILITY_NSI_T:-${NSI_EXTERNAL_NSI_T:-${PFM_EXTERNAL_NSI_T:-10}}}"
NSI_RELIABILITY_QUERY_T="${NSI_RELIABILITY_QUERY_T:-${NSI_EXTERNAL_QUERY_T:-${PFM_EXTERNAL_QUERY_T:-60}}}"

NSI_EXTERNAL_ROOT="${NSI_EXTERNAL_ROOT:-${PFM_EXTERNAL_ROOT:-}}"
NSI_EXTERNAL_ENTRY="${NSI_EXTERNAL_ENTRY:-${PFM_EXTERNAL_ENTRY:-pfm_nsi.cli}}"
NSI_EXTERNAL_OUT_SUBDIR="${NSI_EXTERNAL_OUT_SUBDIR:-${PFM_EXTERNAL_OUT_SUBDIR:-}}"
NSI_EXTERNAL_PREFIX="${NSI_EXTERNAL_PREFIX:-${PFM_EXTERNAL_PREFIX:-pfm_nsi}}"
NSI_EXTERNAL_THRESHOLDS="${NSI_EXTERNAL_THRESHOLDS:-${PFM_EXTERNAL_THRESHOLDS:-0.6,0.7,0.8}}"
NSI_EXTERNAL_MORANS="${NSI_EXTERNAL_MORANS:-${PFM_EXTERNAL_MORANS:-1}}"
NSI_EXTERNAL_SLOPE="${NSI_EXTERNAL_SLOPE:-${PFM_EXTERNAL_SLOPE:-1}}"
NSI_EXTERNAL_RIDGE_LAMBDAS="${NSI_EXTERNAL_RIDGE_LAMBDAS:-${PFM_EXTERNAL_RIDGE_LAMBDAS:-10}}"
NSI_EXTERNAL_STRUCTURES="${NSI_EXTERNAL_STRUCTURES:-${NSI_EXTERNAL_STRUCTURES_CSV:-${PFM_NSI_STRUCTURES_CSV:-}}}"
NSI_EXTERNAL_SPARSE_FRAC="${NSI_EXTERNAL_SPARSE_FRAC:-${PFM_EXTERNAL_SPARSE_FRAC:-}}"
NSI_EXTERNAL_THREADS="${NSI_EXTERNAL_THREADS:-${PFM_EXTERNAL_THREADS:-1}}"
NSI_EXTERNAL_FULLMEM="${NSI_EXTERNAL_FULLMEM:-${PFM_EXTERNAL_FULLMEM:-0}}"
NSI_EXTERNAL_DTYPE="${NSI_EXTERNAL_DTYPE:-${PFM_EXTERNAL_DTYPE:-float32}}"
NSI_EXTERNAL_BLOCK_SIZE="${NSI_EXTERNAL_BLOCK_SIZE:-${PFM_EXTERNAL_BLOCK_SIZE:-2048}}"
NSI_EXTERNAL_KEEP_ALLRHO="${NSI_EXTERNAL_KEEP_ALLRHO:-${PFM_EXTERNAL_KEEP_ALLRHO:-0}}"
NSI_EXTERNAL_KEEP_BETAS="${NSI_EXTERNAL_KEEP_BETAS:-${PFM_EXTERNAL_KEEP_BETAS:-0}}"
NSI_EXTERNAL_KEEP_FC_MAP="${NSI_EXTERNAL_KEEP_FC_MAP:-${PFM_EXTERNAL_KEEP_FC_MAP:-0}}"

if [[ "$NSI_ENABLE" != "1" ]]; then
  echo "[nsi] disabled (NSI_ENABLE=$NSI_ENABLE); skipping"
  exit 0
fi

if [[ "$NSI_USE_EXTERNAL_CLI" != "1" ]]; then
  echo "[nsi] NSI_USE_EXTERNAL_CLI=$NSI_USE_EXTERNAL_CLI; this module currently runs only through the external NSI CLI"
  exit 0
fi
if [[ "$NSI_USABILITY_MODEL" != "1" && "$NSI_RELIABILITY_MODEL" != "1" ]]; then
  echo "[nsi] both NSI_USABILITY_MODEL and NSI_RELIABILITY_MODEL are disabled; skipping"
  exit 0
fi

if [[ -z "$NSI_EXTERNAL_ROOT" ]]; then
  echo "ERROR: NSI_USE_EXTERNAL_CLI=1 but NSI_EXTERNAL_ROOT is empty"
  exit 2
fi

FD_TAG="${NSI_FD_THRESHOLD//./p}"
BASE="${FuncFilePrefix}_${NSI_INPUT_TAG}"
CENSORED_CIFTI="${NSI_INPUT_CIFTI:-${StudyFolder}/${Subject}/func/${FuncDirName}/${NSI_CONCAT_OUT_SUBDIR}/${BASE}_Concatenated+FDlt${FD_TAG}.dtseries.nii}"
SCAN_INFO_JSON="${NSI_SCAN_INFO_JSON:-${StudyFolder}/${Subject}/func/${FuncDirName}/${NSI_CONCAT_OUT_SUBDIR}/ScanInfo.json}"
QC_OUTPUT_DIRNAME="${QC_OUTPUT_DIRNAME:-QC}"
EXTERNAL_OUTDIR_BASE="${NSI_OUTDIR:-${PFM_DATASET_ROOT:-${StudyFolder}/${Subject}/func/${FuncDirName}}/${QC_OUTPUT_DIRNAME}/NSI}"
EXTERNAL_OUTDIR="$EXTERNAL_OUTDIR_BASE"
if [[ -n "${NSI_EXTERNAL_OUT_SUBDIR}" ]]; then
  EXTERNAL_OUTDIR="$EXTERNAL_OUTDIR_BASE/${NSI_EXTERNAL_OUT_SUBDIR}"
fi

if [[ ! -f "$CENSORED_CIFTI" ]]; then
  echo "ERROR: expected censored CIFTI not found for NSI: $CENSORED_CIFTI"
  echo "Run concat stage first (or align NSI_INPUT_TAG/NSI_FD_THRESHOLD with existing outputs)."
  exit 2
fi
# Model duration is the actual input CIFTI duration after frame selection/censoring.
TOTAL_MINUTES="$("$NSI_PYTHON" - "$CENSORED_CIFTI" <<'PYTIME'
import sys
import nibabel as nib
from nibabel.cifti2.cifti2_axes import SeriesAxis
img = nib.load(sys.argv[1])
axis = img.header.get_axis(0)
if not isinstance(axis, SeriesAxis) or axis.unit != 'SECOND' or axis.step <= 0:
    raise ValueError('NSI requires a CIFTI time axis in seconds with positive TR')
if axis.size < 2:
    raise ValueError('NSI requires at least two retained frames')
print('{:.6f}'.format(axis.size * axis.step / 60.0))
PYTIME
)"

mkdir -p "$EXTERNAL_OUTDIR"

EXTERNAL_WORKDIR="$NSI_EXTERNAL_ROOT"
if [[ -d "${NSI_EXTERNAL_ROOT}/pfm-nsi/pfm_nsi" ]]; then
  EXTERNAL_WORKDIR="${NSI_EXTERNAL_ROOT}/pfm-nsi"
fi
if [[ ! -d "${EXTERNAL_WORKDIR}/pfm_nsi" ]]; then
  echo "ERROR: could not locate pfm_nsi package under NSI_EXTERNAL_ROOT='${NSI_EXTERNAL_ROOT}'"
  echo "Checked: ${EXTERNAL_WORKDIR}/pfm_nsi"
  exit 2
fi

echo "[nsi] external NSI CLI root: ${NSI_EXTERNAL_ROOT}"
echo "[nsi] external NSI workdir: ${EXTERNAL_WORKDIR}"
echo "[nsi] external NSI output dir: ${EXTERNAL_OUTDIR}"
echo "[nsi] retained input duration (minutes): ${TOTAL_MINUTES}"

run_external_nsi() {
  local model="$1"
  local nsi_t="$2"
  local query_t="$3"
  local prefix="$4"

  local -a args=(
    -m "$NSI_EXTERNAL_ENTRY" run
    --cifti "$CENSORED_CIFTI"
    --outdir "$EXTERNAL_OUTDIR"
    --prefix "$prefix"
    --nsi-t "$nsi_t"
    --query-t "$query_t"
    --thresholds "$NSI_EXTERNAL_THRESHOLDS"
    --ridge-lambdas "$NSI_EXTERNAL_RIDGE_LAMBDAS"
  )

  if [[ "$model" == "usability" ]]; then
    args+=(--usability)
  elif [[ "$model" == "reliability" ]]; then
    args+=(--reliability)
  else
    echo "ERROR: unknown NSI model '$model'"
    exit 2
  fi

  if [[ "$NSI_EXTERNAL_MORANS" == "1" ]]; then
    args+=(--morans)
  fi
  if [[ "$NSI_EXTERNAL_SLOPE" == "1" ]]; then
    args+=(--slope)
  fi
  if [[ -n "$NSI_EXTERNAL_SPARSE_FRAC" ]]; then
    args+=(--sparse-frac "$NSI_EXTERNAL_SPARSE_FRAC")
  fi
  if [[ -n "$NSI_EXTERNAL_STRUCTURES" ]]; then
    args+=(--structures "$NSI_EXTERNAL_STRUCTURES")
  fi
  args+=(--dtype "$NSI_EXTERNAL_DTYPE" --block-size "$NSI_EXTERNAL_BLOCK_SIZE")
  if [[ "$NSI_EXTERNAL_FULLMEM" == "1" ]]; then
    args+=(--fullmem)
  fi
  if [[ "$NSI_EXTERNAL_KEEP_ALLRHO" == "1" ]]; then
    args+=(--keep-allrho)
  fi
  if [[ "$NSI_EXTERNAL_KEEP_BETAS" == "1" ]]; then
    args+=(--keep-betas)
  fi
  if [[ "$NSI_EXTERNAL_KEEP_FC_MAP" == "1" ]]; then
    args+=(--keep-fc-map)
  fi

  {
    echo "date=$(date --iso-8601=seconds)"
    echo "subject=${Subject}"
    echo "cifti=${CENSORED_CIFTI}"
    echo "python=${NSI_PYTHON}"
    echo "entry=${NSI_EXTERNAL_ENTRY}"
    echo "outdir=${EXTERNAL_OUTDIR}"
    echo "model=${model}"
    echo "prefix=${prefix}"
    echo "nsi_t=${nsi_t}"
    echo "query_t=${query_t}"
    echo "total_concat_minutes=${TOTAL_MINUTES}"
    echo "thresholds=${NSI_EXTERNAL_THRESHOLDS}"
    echo "ridge_lambdas=${NSI_EXTERNAL_RIDGE_LAMBDAS}"
    echo "structures=${NSI_EXTERNAL_STRUCTURES}"
    echo "morans=${NSI_EXTERNAL_MORANS}"
    echo "slope=${NSI_EXTERNAL_SLOPE}"
    echo "sparse_frac=${NSI_EXTERNAL_SPARSE_FRAC}"
    echo "threads=${NSI_EXTERNAL_THREADS}"
    echo "fullmem=${NSI_EXTERNAL_FULLMEM}"
    echo "dtype=${NSI_EXTERNAL_DTYPE}"
    echo "block_size=${NSI_EXTERNAL_BLOCK_SIZE}"
    echo "keep_allrho=${NSI_EXTERNAL_KEEP_ALLRHO}"
    echo "keep_betas=${NSI_EXTERNAL_KEEP_BETAS}"
    echo "keep_fc_map=${NSI_EXTERNAL_KEEP_FC_MAP}"
    printf "argv="
    printf '%q ' "$NSI_PYTHON" "${args[@]}"
    echo
  } > "${EXTERNAL_OUTDIR}/${prefix}_run_config.txt"

  echo "[nsi] running ${model} model: prefix=${prefix} nsi_t=${nsi_t} query_t=${query_t}"
  (
    cd "$EXTERNAL_WORKDIR"
    OMP_NUM_THREADS="$NSI_EXTERNAL_THREADS" \
    OPENBLAS_NUM_THREADS="$NSI_EXTERNAL_THREADS" \
    MKL_NUM_THREADS="$NSI_EXTERNAL_THREADS" \
    NUMEXPR_NUM_THREADS="$NSI_EXTERNAL_THREADS" \
      "$NSI_PYTHON" "${args[@]}"
  )
}

if [[ "$NSI_USABILITY_MODEL" == "1" ]]; then
  # Usability should reflect all available concatenated time.
  FULL_TAG="${TOTAL_MINUTES//./p}"
  run_external_nsi "usability" "$TOTAL_MINUTES" "$TOTAL_MINUTES" "${NSI_EXTERNAL_PREFIX}_usability_full${FULL_TAG}m"
fi

if [[ "$NSI_RELIABILITY_MODEL" == "1" ]]; then
  if awk "BEGIN{exit !($NSI_RELIABILITY_NSI_T < $TOTAL_MINUTES)}"; then
    echo "[nsi] reliability nsi_t (${NSI_RELIABILITY_NSI_T}m) is shorter than full concat (${TOTAL_MINUTES}m); running dedicated reliability pass"
  fi
  REL_TAG="${NSI_RELIABILITY_NSI_T//./p}"
  run_external_nsi "reliability" "$NSI_RELIABILITY_NSI_T" "$NSI_RELIABILITY_QUERY_T" "${NSI_EXTERNAL_PREFIX}_reliability_nsiT${REL_TAG}m"
fi

echo "[nsi] complete"

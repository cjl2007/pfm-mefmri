#!/usr/bin/env bash
# Concat module: concatenate run CIFTIs + FD censoring prep.
set -euo pipefail

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
MEDIR="${3:?missing MEDIR}"
StartSession="${4:?missing StartSession}"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

CONCAT_PYTHON="${CONCAT_PYTHON:-${PFM_PYTHON:-python3}}"
CONCAT_INPUT_TAG="${CONCAT_INPUT_TAG:-${PFM_INPUT_TAG:-OCME+MEICA+MGTR}}"
CONCAT_OUT_SUBDIR="${CONCAT_OUT_SUBDIR:-${PFM_CONCAT_OUT_SUBDIR:-ConcatenatedCiftis}}"
CONCAT_FD_THRESHOLD="${CONCAT_FD_THRESHOLD:-${PFM_FD_THRESHOLD:-0.3}}"
CONCAT_DEMEAN_RUNS="${CONCAT_DEMEAN_RUNS:-${PFM_DEMEAN_RUNS:-1}}"
CONCAT_VAR_NORM_RUNS="${CONCAT_VAR_NORM_RUNS:-${PFM_VAR_NORM_RUNS:-0}}"
CONCAT_VAR_NORM_EPS="${CONCAT_VAR_NORM_EPS:-${PFM_VAR_NORM_EPS:-1e-8}}"
CONCAT_CENSOR_BY_FD="${CONCAT_CENSOR_BY_FD:-${PFM_CENSOR_BY_FD:-1}}"
CONCATENATE_RUNS="${CONCATENATE_RUNS:-${PFM_CONCATENATE_RUNS:-1}}"
CONCAT_SAVE_FD_TXT="${CONCAT_SAVE_FD_TXT:-${PFM_SAVE_FD_TXT:-1}}"
CONCAT_SAVE_SCANIDX_TXT="${CONCAT_SAVE_SCANIDX_TXT:-${PFM_SAVE_SCANIDX_TXT:-1}}"

echo "[concat] start subject=${Subject} input_tag=${CONCAT_INPUT_TAG}"
echo "[concat] output dir: ${StudyFolder}/${Subject}/func/${FuncDirName}/${CONCAT_OUT_SUBDIR}"

"$CONCAT_PYTHON" "$MEDIR/lib/pfm_module.py" \
  --subject "$Subject" \
  --study-folder "$StudyFolder" \
  --mode "prep_only" \
  --start-session "$StartSession" \
  --func-dirname "$FuncDirName" \
  --func-file-prefix "$FuncFilePrefix" \
  --input-tag "$CONCAT_INPUT_TAG" \
  --fd-threshold "$CONCAT_FD_THRESHOLD" \
  --out-subdir "$CONCAT_OUT_SUBDIR" \
  --demean-runs "$CONCAT_DEMEAN_RUNS" \
  --var-norm-runs "$CONCAT_VAR_NORM_RUNS" \
  --var-norm-eps "$CONCAT_VAR_NORM_EPS" \
  --censor-by-fd "$CONCAT_CENSOR_BY_FD" \
  --concatenate-runs "$CONCATENATE_RUNS" \
  --save-fd-txt "$CONCAT_SAVE_FD_TXT" \
  --save-scanidx-txt "$CONCAT_SAVE_SCANIDX_TXT"

echo "[concat] complete"

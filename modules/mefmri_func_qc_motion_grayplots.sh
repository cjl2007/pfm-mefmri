#!/usr/bin/env bash
# Generate Motion & GrayPlots QC reports from run-level CIFTI stages.

set -euo pipefail

if [[ "$#" -lt 4 || "$#" -gt 6 ]]; then
  echo "Usage: mefmri_func_qc_motion_grayplots.sh <MEDIR> <Subject> <StudyFolder> <StartSession> [FuncDirName] [FuncFilePrefix]" >&2
  exit 2
fi

MEDIR="$1"
Subject="$2"
StudyFolder="$3"
StartSession="$4"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

[[ "${MEDIR:0:1}" == "/" ]] || MEDIR="$(cd "$MEDIR" && pwd)"
[[ "${StudyFolder:0:1}" == "/" ]] || StudyFolder="$(cd "$StudyFolder" && pwd)"

Subdir="$StudyFolder/$Subject"
ReportPy="${MOTION_GRAYPLOTS_QC_SCRIPT:-$MEDIR/lib/mefmri_qc_motion_grayplots.py}"
PythonBin="${MOTION_GRAYPLOTS_QC_PYTHON:-${PIPELINE_PYTHON:-python3}}"
Enable="${MOTION_GRAYPLOTS_QC_ENABLE:-1}"
Strict="${MOTION_GRAYPLOTS_QC_STRICT:-0}"
MaxRows="${MOTION_GRAYPLOTS_QC_MAX_ROWS:-1200}"
QCFCSample="${MOTION_GRAYPLOTS_QC_QCFC_SAMPLE_SIZE:-300}"
FDThreshold="${MOTION_GRAYPLOTS_QC_FD_THRESHOLD:-0.3}"

if [[ "$Enable" != "1" ]]; then
  echo "[qc-motion-grayplots] skipped (MOTION_GRAYPLOTS_QC_ENABLE=$Enable)"
  exit 0
fi

if [[ ! -f "$ReportPy" ]]; then
  echo "[qc-motion-grayplots] WARNING: missing report script: $ReportPy" >&2
  [[ "$Strict" == "1" ]] && exit 2 || exit 0
fi

echo "[qc-motion-grayplots] generating reports for func/${FuncDirName}"
set +e
"$PythonBin" "$ReportPy" \
  --subject-dir "$Subdir" \
  --func-name "$FuncDirName" \
  --func-prefix "$FuncFilePrefix" \
  --start-session "$StartSession" \
  --max-grayplot-rows "$MaxRows" \
  --qcfc-sample-size "$QCFCSample" \
  --fd-threshold "$FDThreshold"
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
  echo "[qc-motion-grayplots] WARNING: report generation failed with exit ${rc}" >&2
  [[ "$Strict" == "1" ]] && exit "$rc" || exit 0
fi

echo "[qc-motion-grayplots] complete"

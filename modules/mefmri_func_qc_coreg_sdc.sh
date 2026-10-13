#!/usr/bin/env bash
# Generate fMRIPrep-style functional SDC/coregistration QC reports.

set -euo pipefail

if [[ "$#" -lt 4 || "$#" -gt 8 ]]; then
  echo "Usage: mefmri_func_qc_coreg_sdc.sh <MEDIR> <Subject> <StudyFolder> <StartSession> [FuncDirName] [FuncFilePrefix] [FuncXfmsDir] [DistortionMode]" >&2
  exit 2
fi

MEDIR="$1"
Subject="$2"
StudyFolder="$3"
StartSession="$4"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"
FuncXfmsDir="${7:-${FUNC_XFMS_DIRNAME:-$FuncDirName}}"
DistortionMode="${8:-${DISTORTION_CORRECTION_MODE:-auto}}"

[[ "${MEDIR:0:1}" == "/" ]] || MEDIR="$(cd "$MEDIR" && pwd)"
[[ "${StudyFolder:0:1}" == "/" ]] || StudyFolder="$(cd "$StudyFolder" && pwd)"

Subdir="$StudyFolder/$Subject"
ReportPy="${COREG_SDC_QC_SCRIPT:-$MEDIR/lib/mefmri_qc_coreg_sdc_report.py}"
PythonBin="${COREG_SDC_QC_PYTHON:-${PIPELINE_PYTHON:-python3}}"
Enable="${COREG_SDC_QC_ENABLE:-1}"
Strict="${COREG_SDC_QC_STRICT:-0}"

if [[ "$Enable" != "1" ]]; then
  echo "[qc-coreg-sdc] skipped (COREG_SDC_QC_ENABLE=$Enable)"
  exit 0
fi

if [[ ! -f "$ReportPy" ]]; then
  echo "[qc-coreg-sdc] WARNING: missing report script: $ReportPy" >&2
  [[ "$Strict" == "1" ]] && exit 2 || exit 0
fi

echo "[qc-coreg-sdc] generating reports for func/${FuncDirName}, xfms/${FuncXfmsDir}, mode=${DistortionMode}"
set +e
"$PythonBin" "$ReportPy" \
  --subject-dir "$Subdir" \
  --func-name "$FuncDirName" \
  --xfms-name "$FuncXfmsDir" \
  --func-prefix "$FuncFilePrefix" \
  --distortion-mode "$DistortionMode"
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
  echo "[qc-coreg-sdc] WARNING: report generation failed with exit ${rc}" >&2
  [[ "$Strict" == "1" ]] && exit "$rc" || exit 0
fi

echo "[qc-coreg-sdc] complete"

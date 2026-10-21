#!/usr/bin/env bash
# Optional post-PFM TANS entrypoint. Scientific logic remains in pfm-tans.
set -euo pipefail

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
_MEDIR="${3:?missing MEDIR}"
_StartSession="${4:?missing StartSession}"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
_FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

SubjectDir="${StudyFolder}/${Subject}"
TANS_CODE_ROOT="${TANS_CODE_ROOT:?TANS_CODE_ROOT is required}"
TANS_CONFIG_FILE="${TANS_CONFIG_FILE:?TANS_CONFIG_FILE is required when TANS_ENABLE=1}"
TANS_PYTHON="${TANS_PYTHON:-python3}"
TANS_OUTPUT_DIRNAME="${TANS_OUTPUT_DIRNAME:-tans}"
if [[ -n "${TANS_OUTPUT_DIR:-}" ]]; then
  echo "ERROR: TANS_OUTPUT_DIR is no longer supported. Use TANS_OUTPUT_DIRNAME (a folder name under SubjectDir, default: tans)." >&2
  exit 2
fi
case "$TANS_OUTPUT_DIRNAME" in
  .|..|*/*|*\\*|*$'\n'*|*$'\r'*)
    echo "ERROR: TANS_OUTPUT_DIRNAME must be a single folder name under SubjectDir, not a path." >&2
    exit 2
    ;;
esac
TansOutputDir="${SubjectDir}/${TANS_OUTPUT_DIRNAME}"
if [[ -z "${TANS_PFM_DIR:-}" ]]; then
  TANS_PFM_DIR="${SubjectDir}/func/${FuncDirName}/PFM/RidgePFM"
  if [[ ! -d "$TANS_PFM_DIR" && -d "${SubjectDir}/func/${FuncDirName}/PFM/RidgeFusion" ]]; then
    TANS_PFM_DIR="${SubjectDir}/func/${FuncDirName}/PFM/RidgeFusion"
    echo "[tans] using existing RidgeFusion output; YAML input filenames must match that output."
  fi
fi
TANS_WORKERS="${TANS_WORKERS:-6}"
TANS_ROBUSTNESS_ENABLE="${TANS_ROBUSTNESS_ENABLE:-1}"
TANS_TARGET_NETWORK="${TANS_TARGET_NETWORK:-}"
TANS_COARSE_GRID_DENSITY="${TANS_COARSE_GRID_DENSITY:-0.85}"
TANS_TARGET_SEARCH_SPACE="${TANS_TARGET_SEARCH_SPACE:-glasser_predefined_plus_subject_specific_mfg}"
TANS_INCLUDE_SM_HAND_CALIBRATION="${TANS_INCLUDE_SM_HAND_CALIBRATION:-1}"
TANS_SM_HAND_COIL_MODEL="${TANS_SM_HAND_COIL_MODEL:-}"
TANS_TREATMENT_COIL_MODEL="${TANS_TREATMENT_COIL_MODEL:-}"
TANS_CALIBRATION_REGISTRY="${TANS_CALIBRATION_REGISTRY:-}"
TANS_SM_HAND_DIDT_CALIBRATION="${TANS_SM_HAND_DIDT_CALIBRATION:-}"
TANS_TREATMENT_DIDT_CALIBRATION="${TANS_TREATMENT_DIDT_CALIBRATION:-}"
TANS_ALLOW_SIMILAR_COIL_CALIBRATION="${TANS_ALLOW_SIMILAR_COIL_CALIBRATION:-0}"
TANS_CALIBRATION_SUBSTITUTION_JUSTIFICATION="${TANS_CALIBRATION_SUBSTITUTION_JUSTIFICATION:-}"
TANS_MAXIMUM_TREATMENT_RMT_RATIO="${TANS_MAXIMUM_TREATMENT_RMT_RATIO-1.20}"
TANS_HANDLE_ORIENTATION_PREFERENCE="${TANS_HANDLE_ORIENTATION_PREFERENCE:-posterior_inferior}"
TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT="${TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT:-1}"
TANS_START_FROM="${TANS_START_FROM:-headmodel}"
TANS_STOP_AFTER="${TANS_STOP_AFTER:-reports}"
TANS_RESUME="${TANS_RESUME:-1}"

[[ -f "${SubjectDir}/anat/T1w/T1w.nii.gz" ]] || {
  echo "ERROR: missing native TANS T1w: ${SubjectDir}/anat/T1w/T1w.nii.gz" >&2
  exit 2
}
TansNativeT2w="${SubjectDir}/anat/T2w/T2w.nii.gz"
if [[ ! -f "$TansNativeT2w" ]]; then
  echo "[tans] WARNING: native T2w unavailable; falling back to a T1w-only CHARM head model" >&2
fi
[[ -f "$TANS_CONFIG_FILE" ]] || {
  echo "ERROR: missing TANS config: $TANS_CONFIG_FILE" >&2
  exit 2
}
[[ -f "$TANS_CODE_ROOT/run_pipeline.py" ]] || {
  echo "ERROR: missing TANS Python driver: $TANS_CODE_ROOT/run_pipeline.py" >&2
  exit 2
}
[[ -d "$TANS_PFM_DIR" ]] || {
  echo "ERROR: missing PFM output: $TANS_PFM_DIR" >&2
  exit 2
}
[[ "$TANS_WORKERS" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: TANS_WORKERS must be positive" >&2
  exit 2
}
[[ "$TANS_ROBUSTNESS_ENABLE" =~ ^[01]$ ]] || {
  echo "ERROR: TANS_ROBUSTNESS_ENABLE must be 0 or 1" >&2
  exit 1
}
[[ "$TANS_INCLUDE_SM_HAND_CALIBRATION" =~ ^[01]$ ]] || {
  echo "ERROR: TANS_INCLUDE_SM_HAND_CALIBRATION must be 0 or 1" >&2
  exit 2
}
[[ "$TANS_ALLOW_SIMILAR_COIL_CALIBRATION" =~ ^[01]$ ]] || {
  echo "ERROR: TANS_ALLOW_SIMILAR_COIL_CALIBRATION must be 0 or 1" >&2
  exit 2
}
[[ "$TANS_HANDLE_ORIENTATION_PREFERENCE" =~ ^(as_simulated|posterior_inferior)$ ]] || {
  echo "ERROR: TANS_HANDLE_ORIENTATION_PREFERENCE must be as_simulated or posterior_inferior" >&2
  exit 2
}
[[ "$TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT" =~ ^[01]$ ]] || {
  echo "ERROR: TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT must be 0 or 1" >&2
  exit 2
}

command=(
  "$TANS_PYTHON" "$TANS_CODE_ROOT/run_pipeline.py" "$TANS_CONFIG_FILE"
  --subject-dir "$SubjectDir"
  --func-dir "$FuncDirName"
  --pfm-dir "$TANS_PFM_DIR"
  --output "$TansOutputDir"
  --workers "$TANS_WORKERS"
  --robustness-enabled "$TANS_ROBUSTNESS_ENABLE"
  --coarse-grid-density "$TANS_COARSE_GRID_DENSITY"
  --handle-orientation-preference "$TANS_HANDLE_ORIENTATION_PREFERENCE"
  --coil-180-degree-magnitude-equivalent "$TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT"
  --start-from "$TANS_START_FROM"
  --stop-after "$TANS_STOP_AFTER"
)
# Unspecified hardware stays with the study YAML; do not inject a site-specific
# installation or a calibration for a different coil into an otherwise valid config.
if [[ -n "$TANS_SM_HAND_COIL_MODEL" ]]; then
  command+=(--sm-hand-coil "$TANS_SM_HAND_COIL_MODEL")
fi
if [[ -n "$TANS_TREATMENT_COIL_MODEL" ]]; then
  command+=(--treatment-coil "$TANS_TREATMENT_COIL_MODEL")
fi
if [[ -n "$TANS_MAXIMUM_TREATMENT_RMT_RATIO" ]]; then
  command+=(--maximum-treatment-rmt-ratio "$TANS_MAXIMUM_TREATMENT_RMT_RATIO")
else
  command+=(--disable-treatment-rmt-limit)
fi
if [[ -n "$TANS_TARGET_NETWORK" ]]; then
  command+=(--target-network "$TANS_TARGET_NETWORK")
fi
if [[ -n "$TANS_TARGET_SEARCH_SPACE" ]]; then
  command+=(--target-search-space "$TANS_TARGET_SEARCH_SPACE")
fi
if [[ -n "$TANS_CALIBRATION_REGISTRY" ]]; then
  command+=(--calibration-registry "$TANS_CALIBRATION_REGISTRY")
fi
if [[ -n "$TANS_SM_HAND_DIDT_CALIBRATION" ]]; then
  command+=(--sm-hand-didt-calibration "$TANS_SM_HAND_DIDT_CALIBRATION")
fi
if [[ -n "$TANS_TREATMENT_DIDT_CALIBRATION" ]]; then
  command+=(--treatment-didt-calibration "$TANS_TREATMENT_DIDT_CALIBRATION")
fi
if [[ "$TANS_ALLOW_SIMILAR_COIL_CALIBRATION" == "1" ]]; then
  command+=(--allow-similar-coil-calibration)
fi
if [[ -n "$TANS_CALIBRATION_SUBSTITUTION_JUSTIFICATION" ]]; then
  command+=(--calibration-substitution-justification \
    "$TANS_CALIBRATION_SUBSTITUTION_JUSTIFICATION")
fi
if [[ "$TANS_INCLUDE_SM_HAND_CALIBRATION" == "0" ]]; then
  command+=(--skip-sm-hand-calibration)
fi
if [[ "$TANS_RESUME" == "1" ]]; then
  command+=(--resume)
fi

echo "[tans] native T1w=${SubjectDir}/anat/T1w/T1w.nii.gz"
if [[ -f "$TansNativeT2w" ]]; then
  echo "[tans] native T2w=${TansNativeT2w}"
else
  echo "[tans] native T2w=unavailable (T1w-only fallback)"
fi
echo "[tans] PFM input root=${TANS_PFM_DIR}"
echo "[tans] output=${TansOutputDir}"
echo "[tans] treatment network=${TANS_TARGET_NETWORK:-none; calibration only}"
echo "[tans] coarse grid density=${TANS_COARSE_GRID_DENSITY}"
echo "[tans] placement-error robustness=${TANS_ROBUSTNESS_ENABLE} (0: nominal fine winner, not evaluated)"
echo "[tans] target search space=${TANS_TARGET_SEARCH_SPACE:-default SM-Hand precentral zone}"
echo "[tans] include SM-Hand calibration=${TANS_INCLUDE_SM_HAND_CALIBRATION}"
echo "[tans] SM-Hand coil model=${TANS_SM_HAND_COIL_MODEL:-study YAML}"
echo "[tans] treatment coil model=${TANS_TREATMENT_COIL_MODEL:-study YAML}"
echo "[tans] calibration registry=${TANS_CALIBRATION_REGISTRY:-pfm-tans default}"
echo "[tans] SM-Hand dI/dt calibration=${TANS_SM_HAND_DIDT_CALIBRATION:-study YAML or exact-match resolution}"
echo "[tans] treatment dI/dt calibration=${TANS_TREATMENT_DIDT_CALIBRATION:-study YAML or exact-match resolution; ambiguity stops for user selection}"
echo "[tans] allow similar-coil calibration=${TANS_ALLOW_SIMILAR_COIL_CALIBRATION}"
echo "[tans] maximum treatment/rMT MSO ratio=${TANS_MAXIMUM_TREATMENT_RMT_RATIO:-disabled}"
echo "[tans] handle orientation preference=${TANS_HANDLE_ORIENTATION_PREFERENCE}"
echo "[tans] coil 180-degree |E| equivalence asserted=${TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT}"
printf '[tans] command:'
printf ' %q' "${command[@]}"
printf '\n'
"${command[@]}"

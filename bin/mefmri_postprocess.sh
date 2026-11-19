#!/usr/bin/env bash
# Independent optional consumers. Never creates or prepares a PFM dataset.
set -euo pipefail
SubjectDir="$1"; ConfigFile="$2"
set -a
source "$ConfigFile"
set +a
: "${TANS_ENABLE:=0}" "${TASK_MODEL_ENABLE:=0}"
for flag in "$TANS_ENABLE" "$TASK_MODEL_ENABLE"; do
  [[ "$flag" == 0 || "$flag" == 1 ]] || { echo 'ERROR: optional enable switches must be 0 or 1' >&2; exit 2; }
done
[[ "$TANS_ENABLE$TASK_MODEL_ENABLE" != 00 ]] || exit 0
if [[ "$TASK_MODEL_ENABLE" == 1 ]]; then
  [[ -d "${TASK_MODEL_CONFIG_DIR:-}" ]] || { echo 'ERROR: TASK_MODEL_CONFIG_DIR is missing' >&2; exit 2; }
  declare -p TASK_MODEL_TASKS >/dev/null 2>&1 || { echo 'ERROR: TASK_MODEL_TASKS must be specified' >&2; exit 2; }
  requested_tasks=("${TASK_MODEL_TASKS[@]}")
  [[ "${#requested_tasks[@]}" -gt 0 ]] || { echo 'ERROR: TASK_MODEL_TASKS must not be empty' >&2; exit 2; }
  declare -A task_seen=()
  for task in "${requested_tasks[@]}"; do
    [[ "$task" =~ ^[A-Za-z0-9_-]+$ && -z "${task_seen[$task]:-}" ]] || { echo "ERROR: invalid or duplicate task name: $task" >&2; exit 2; }
    task_seen[$task]=1
    [[ -f "$TASK_MODEL_CONFIG_DIR/$task.yaml" ]] || { echo "ERROR: missing task recipe: $TASK_MODEL_CONFIG_DIR/$task.yaml" >&2; exit 2; }
  done
  export TASK_MODEL_TASKS_JSON="$("${TASK_MODEL_PYTHON:-${PIPELINE_PYTHON:-python3}}" -c 'import json,sys;print(json.dumps(sys.argv[1:]))' "${requested_tasks[@]}")"
fi
source "$EnvironmentScript"
Subject="$(basename "$SubjectDir")"; StudyFolder="$(dirname "$SubjectDir")"
logdir="$SubjectDir/QC/Postprocessing/ModuleLogs"; mkdir -p "$logdir"
stamp="$(date +%Y%m%d_%H%M%S)_$$"
run_optional() {
  local label="$1" module="$2"
  [[ -f "$module" ]] || { echo "ERROR: missing $label module: $module" >&2; return 2; }
  echo "[postprocessing] $label log: $logdir/${stamp}_${label}.log"
  bash "$module" "$Subject" "$StudyFolder" "$MEDIR" "${START_SESSION:-1}" "${FUNC_DIRNAME[0]:-rest}" "${FUNC_FILE_PREFIX[0]:-Rest}" > "$logdir/${stamp}_${label}.log" 2>&1 || return $?
  echo "[postprocessing] $label complete"
}
# The consumers do not depend on each other's outputs.
status=0
if [[ "$TANS_ENABLE" == 1 ]]; then
  run_optional tans "${FUNC_TANS_MODULE:-$MEDIR/modules/mefmri_func_tans.sh}" || status=$?
fi
if [[ "$TASK_MODEL_ENABLE" == 1 ]]; then
  run_optional task_model "${FUNC_TASK_MODEL_MODULE:-$MEDIR/modules/mefmri_func_task_model.sh}" || status=$?
fi
exit "$status"

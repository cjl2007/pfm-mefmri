#!/usr/bin/env bash
set -euo pipefail
echo "[task_model] Experimental — work in progress"
Subject="$1"; StudyFolder="$2"; MEDIR="$3"; StartSession="$4"
exec "${TASK_MODEL_PYTHON:-${PIPELINE_PYTHON:-python3}}" "$MEDIR/lib/mefmri_task_model.py" \
  --subject-dir "$StudyFolder/$Subject" --registry "${TASK_MODEL_CONFIG_DIR:-$MEDIR/config/task_models}" \
  --start-session "$StartSession" --mode fit

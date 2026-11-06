#!/usr/bin/env bash
# Cheap raw-session feasibility check; retained-frame eligibility is checked at PFM entry.
[[ "${PFM_ENABLE:-0}" == 1 && "${PFM_RF_CONFIDENCE_ENABLE:-0}" == 1 ]] || return 0
case "${START_FROM_MODULE:-}" in tans|pfm_update|infomap_update) return 0 ;; esac
case "${PFM_STRATEGY:-ridgepfm}" in
  ridgepfm|RidgePFM|ridge_fusion|RidgeFusion|ridge-fusion|all|All|ALL) ;;
  *) return 0 ;;
esac
case "${STOP_AFTER_MODULE:-}" in
  validate|pfm|pfm_update|infomap_update|tans|"") ;;
  *) return 0 ;;
esac
# An explicitly supplied CIFTI can have a different session selection.
[[ -z "${PFM_INPUT_CIFTI:-}" || -n "${PFM_INPUT_TASK_SPECS_JSON:-}" ]] || return 0
pfm_confidence_raw_preflight() {
  local subject_dir="$1" func_dir="$2" start_session="$3" session count=0 number
  if [[ -n "${PFM_INPUT_TASK_SPECS_JSON:-}" ]]; then
    count="$("${PFM_PYTHON:-${PIPELINE_PYTHON:-python3}}" - "$subject_dir" "$start_session" "$PFM_INPUT_TASK_SPECS_JSON" <<'PY'
import json
from pathlib import Path
import re
import sys
subject, start, specs = Path(sys.argv[1]), int(sys.argv[2]), json.loads(sys.argv[3])
sessions = set()
for spec in specs:
    for root in (subject / 'func/unprocessed' / spec['task'], subject / 'func' / spec['task']):
        for directory in root.glob('session_*'):
            match = re.fullmatch(r'session_([0-9]+)', directory.name)
            if match and int(match.group(1)) >= start and any(p.is_dir() for p in directory.glob('run_*')):
                sessions.add(int(match.group(1)))
print(len(sessions))
PY
)" || return 2
  else
  for session in "$subject_dir/func/unprocessed/$func_dir"/session_*; do
    [[ -d "$session" && "$session" =~ /session_([0-9]+)$ ]] || continue
    number="${BASH_REMATCH[1]}"
    (( 10#$number >= 10#$start_session )) || continue
    local run
    for run in "$session"/run_*; do
      [[ -d "$run" && "$run" =~ /run_[0-9]+$ ]] || continue
      count=$((count + 1))
      break
    done
  done
  fi
  if (( count < 4 )); then
    echo "ERROR: PFM_RF_CONFIDENCE_ENABLE=1 requires at least four eligible sessions (two per half); found $count raw sessions with runs from session $start_session in func/$func_dir. Multiple runs in one session do not count separately. Set PFM_RF_CONFIDENCE_ENABLE=0 to run RidgePFM without confidence analysis." >&2
    return 2
  fi
  echo "[validate] RidgePFM confidence: $count raw sessions; at least four must retain ${PFM_RF_REPEAT_MIN_FRAMES:-30} frames each after censoring (checked at PFM entry)."
}
pfm_confidence_raw_preflight "$1" "$2" "$3"

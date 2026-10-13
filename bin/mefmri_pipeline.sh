#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: $0 <SubjectDir> [ConfigFile]"
  echo "Example: $0 /path/to/study/ME01"
  exit 2
fi

SubjectDir="$1"
ConfigFileArg="${2:-}"

if [ "${SubjectDir: -1}" = "/" ]; then
  SubjectDir="${SubjectDir%?}"
fi
if [ ! -d "$SubjectDir" ]; then
  echo "ERROR: subject directory does not exist: $SubjectDir"
  exit 2
fi

Subject="$(basename "$SubjectDir")"
StudyFolder="$(dirname "$SubjectDir")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Give copied study configs the installation root before they define resource paths.
PFM_MEFMRI_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEFAULT_CONFIG="$PFM_MEFMRI_ROOT/config/mefmri_wrapper_config.sh"
ConfigFile="${ConfigFileArg:-${CONFIG_FILE:-$DEFAULT_CONFIG}}"

if [ ! -f "$ConfigFile" ]; then
  echo "ERROR: missing config file: $ConfigFile"
  exit 2
fi

set -a
source "$ConfigFile"
set +a

if [[ -n "${PIPELINE_TASK_OVERRIDE:-}" ]]; then
  FUNC_DIRNAME="$PIPELINE_TASK_OVERRIDE"
fi
if [[ -n "${PIPELINE_TASK_XFMS_OVERRIDE:-}" ]]; then
  FUNC_XFMS_DIRNAME="$PIPELINE_TASK_XFMS_OVERRIDE"
fi

# Optional config-specified FreeSurfer license path.
if [[ -n "${FS_LICENSE_FILE:-}" ]]; then
  export FS_LICENSE="$FS_LICENSE_FILE"
fi

# Core paths and wrapper controls.
: "${MEDIR:=$(cd "$SCRIPT_DIR/.." && pwd)}"
: "${EnvironmentScript:=$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh}"
: "${START_SESSION:=1}"
: "${START_FROM_MODULE:=validate}" # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|denoise|mgtr|vol2surf|concat|nsi|pfm|pfm_update|tans
: "${STOP_AFTER_MODULE:=}" # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|denoise|mgtr|vol2surf|concat|nsi|pfm|pfm_update|tans
: "${PIPELINE_SKIP_ANAT:=0}" # 0|1 (used by FUNC_DIRNAME=All parent loop to avoid rerunning anatomy per task)

# Optional post-config overrides (used internally by the All-tasks driver).
: "${PIPELINE_START_FROM_MODULE_OVERRIDE:=}"
: "${PIPELINE_STOP_AFTER_MODULE_OVERRIDE:=}"

# Module entrypoints.
: "${VALIDATE_MODULE:=$MEDIR/modules/mefmri_validate_inputs.sh}"
: "${ANAT_HCP_MODULE:=$MEDIR/modules/mefmri_anat_hcp.sh}"
: "${ANAT_CHARM_MODULE:=$MEDIR/modules/mefmri_anat_charm.sh}"
: "${FUNC_FIELDMAPS_MODULE:=$MEDIR/modules/mefmri_func_fieldmaps.sh}"
: "${FUNC_MEDIC_MODULE:=$MEDIR/modules/mefmri_func_medic.sh}"
: "${FUNC_MEDIC_COREG_MODULE:=$MEDIR/modules/mefmri_func_coreg_medic.sh}"
: "${FUNC_MEDIC_HEADMOTION_MODULE:=$MEDIR/modules/mefmri_func_headmotion_medic.sh}"
: "${FUNC_QC_COREG_SDC_MODULE:=$MEDIR/modules/mefmri_func_qc_coreg_sdc.sh}"
: "${FUNC_QC_MOTION_GRAYPLOTS_MODULE:=$MEDIR/modules/mefmri_func_qc_motion_grayplots.sh}"
: "${FUNC_COREG_MODULE:=$MEDIR/modules/mefmri_func_coreg.sh}"
: "${FUNC_HEADMOTION_MODULE:=$MEDIR/modules/mefmri_func_headmotion.sh}"
: "${FUNC_MEICA_MODULE:=$MEDIR/modules/mefmri_func_meica.sh}"
: "${FUNC_SINGLEECHO_MODULE:=$MEDIR/modules/mefmri_func_singleecho.sh}"
: "${FUNC_AROMA_MODULE:=$MEDIR/modules/mefmri_func_aroma.sh}"
: "${FUNC_MGTR_MODULE:=$MEDIR/modules/mefmri_func_mgtr.sh}"
: "${FUNC_VOL2SURF_MODULE:=$MEDIR/modules/mefmri_func_vol2surf.sh}"
: "${FUNC_CONCAT_MODULE:=$MEDIR/modules/mefmri_func_concat.sh}"
: "${FUNC_NSI_MODULE:=$MEDIR/modules/mefmri_func_nsi.sh}"
: "${FUNC_PFM_MODULE:=$MEDIR/modules/mefmri_func_pfm.sh}"
: "${FUNC_PFM_INFOMAP_UPDATE_MODULE:=$MEDIR/modules/mefmri_func_pfm_infomap_update.sh}"
: "${FUNC_CRAWLING_SEED_FC_MODULE:=$MEDIR/modules/mefmri_func_crawling_seed_fc.sh}"
: "${FUNC_TANS_MODULE:=$MEDIR/modules/mefmri_func_tans.sh}"

# Global processing knobs.
: "${CHARM_BIN:=}" # optional explicit path passed to CHARM module
: "${DOF:=6}"
: "${AtlasTemplate:=$MEDIR/res0urces/MNI152_T1_2mm.nii.gz}"
: "${AtlasSpace:=T1w}"
: "${MEPCA:=kundu}"
: "${MaxIterations:=500}"
: "${MaxRestarts:=5}"
: "${CONCAT_ENABLE:=1}" # 0|1
: "${NSI_ENABLE:=1}" # 0|1
: "${PFM_ENABLE:=1}" # 0|1
: "${PFM_INFOMAP_UPDATE_ENABLE:=0}" # 0|1
: "${TANS_ENABLE:=0}" # 0|1
: "${CRAWLING_SEED_FC_ENABLE:=1}" # 0|1
: "${RUN_CONFIG_SNAPSHOT:=1}" # 0|1
: "${QC_OUTPUT_DIRNAME:=QC}"
: "${QC_VALIDATION_ENABLE:=1}" # 0|1
: "${FUNC_NOFIELDMAP_MODE:=0}" # 0|1
: "${DISTORTION_CORRECTION_MODE:=topup}" # topup|direct_b0|medic|none
: "${COREG_SDC_QC_ENABLE:=1}" # 0|1
: "${COREG_SDC_QC_STRICT:=0}" # 0|1
: "${MOTION_GRAYPLOTS_QC_ENABLE:=1}" # 0|1
: "${MOTION_GRAYPLOTS_QC_STRICT:=0}" # 0|1
: "${MOTION_GRAYPLOTS_QC_FD_THRESHOLD:=0.3}"
: "${PROCESSING_MODE:=auto}" # auto|multi_echo|single_echo
: "${MULTI_ECHO_DENOISE_METHOD:=meica}" # meica|aroma|acompcor
: "${SINGLE_ECHO_DENOISE_METHOD:=acompcor}" # aroma|acompcor
: "${SINGLE_ECHO_ECHO_INDEX:=1}"
: "${AROMA_NSI_RECLASSIFY_ENABLE:=0}" # 0|1; future experimental path, classic AROMA by default
: "${AROMA_NSI_THRESHOLD:=0.05}"

case "$DISTORTION_CORRECTION_MODE" in
  topup|direct_b0|fieldmap|gre_fieldmap) ;;
  none)
    FUNC_NOFIELDMAP_MODE=1
    ;;
  medic)
    # MEDIC requires per-echo phase companions and provides its own
    # coreg/headmotion modules after displacement-map generation.
    ;;
  *)
    echo "ERROR: DISTORTION_CORRECTION_MODE must be topup, direct_b0, medic, or none (got '$DISTORTION_CORRECTION_MODE')" >&2
    exit 2
    ;;
esac

# Functional naming/outputs.
: "${VOL2SURF_INPUTS:=}"
: "${FUNC_DIRNAME:=rest}"
: "${FUNC_FILE_PREFIX:=Rest}"
: "${PIPELINE_TASK_LOOP_CHILD:=0}"
if [[ -z "${FUNC_XFMS_DIRNAME:-}" ]]; then
  FUNC_XFMS_DIRNAME="$FUNC_DIRNAME"
fi

has_topup_fieldmaps() {
  local fm_rel="${FM_RAW_DIR_REL:-func/unprocessed/${FUNC_DIRNAME}/field_maps}"
  local fm_dir="$SubjectDir/$fm_rel"
  compgen -G "$fm_dir/AP_S*_R*.nii.gz" >/dev/null || return 1
  compgen -G "$fm_dir/PA_S*_R*.nii.gz" >/dev/null || return 1
}

medic_phase_inputs_complete() {
  local root="$SubjectDir/func/unprocessed/$FUNC_DIRNAME"
  local session_dir run_dir s r mag phase
  local -a mags=()
  local found_run=0
  [[ -d "$root" ]] || return 1

  while IFS= read -r session_dir; do
    s="${session_dir##*/session_}"
    [[ "$s" =~ ^[0-9]+$ ]] || continue
    (( s >= START_SESSION )) || continue
    while IFS= read -r run_dir; do
      r="${run_dir##*/run_}"
      [[ "$r" =~ ^[0-9]+$ ]] || continue
      mapfile -t mags < <(find -L "$run_dir" -maxdepth 1 -type f -name "${FUNC_FILE_PREFIX}_S${s}_R${r}_E*.nii.gz" ! -name '*_phase.nii.gz' | sort -V)
      [[ "${#mags[@]}" -ge 2 ]] || return 1
      for mag in "${mags[@]}"; do
        phase="${mag%.nii.gz}_phase.nii.gz"
        [[ -f "$phase" ]] || return 1
      done
      found_run=1
    done < <(find -L "$session_dir" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
  done < <(find -L "$root" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V)

  [[ "$found_run" == "1" ]]
}

has_reusable_topup_outputs() {
  local fm_rel="${FM_OUT_DIR_REL:-func/${FUNC_DIRNAME}/field_maps}"
  local fm_dir="$SubjectDir/$fm_rel"
  local output mode
  [[ -r "$fm_dir/FIELDMAP_MODE.txt" ]] || return 1
  mode="$(cat "$fm_dir/FIELDMAP_MODE.txt")"
  # The fieldmaps module records TOPUP as "measured", direct B0 as
  # "direct_b0", and no-fieldmap placeholders as "none".
  [[ "$mode" == "measured" || "$mode" == "topup" ]] || return 1
  for output in Avg_FM_rads_acpc.nii.gz Avg_FM_mag_acpc.nii.gz \
                Avg_FM_mag_acpc_brain.nii.gz acqparams.txt; do
    [[ -s "$fm_dir/$output" ]] || return 1
  done
}

resolve_medic_resume_mode() {
  local root="$SubjectDir/func/$FUNC_DIRNAME"
  local session_dir run_dir s r dmap first_missing=""
  local total=0 present=0 complete=0
  # Match the input run inventory of the stage being resumed. Coreg can
  # recreate processed run directories; headmotion consumes those directories.
  if [[ "$START_FROM_MODULE" == "coreg" && -d "$SubjectDir/func/unprocessed/$FUNC_DIRNAME" ]]; then
    root="$SubjectDir/func/unprocessed/$FUNC_DIRNAME"
  fi
  if [[ -d "$root" ]]; then
    while IFS= read -r session_dir; do
      s="${session_dir##*/session_}"
      [[ "$s" =~ ^[0-9]+$ ]] || continue
      (( s >= START_SESSION )) || continue
      while IFS= read -r run_dir; do
        r="${run_dir##*/run_}"
        [[ "$r" =~ ^[0-9]+$ ]] || continue
        total=$((total + 1))
        dmap="$SubjectDir/func/$FUNC_DIRNAME/session_${s}/run_${r}/MEDIC/${FUNC_FILE_PREFIX}_S${s}_R${r}_displacementmaps.nii"
        if [[ -e "$dmap" || -L "$dmap" || -e "${dmap}.gz" || -L "${dmap}.gz" ]]; then
          present=$((present + 1))
        fi
        if [[ -s "$dmap" ]]; then
          complete=$((complete + 1))
        elif [[ -z "$first_missing" ]]; then
          first_missing="$dmap"
        fi
      done < <(find -L "$session_dir" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
    done < <(find -L "$root" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V)
  fi

  if (( total > 0 && complete == total )); then
    echo "Resuming with existing MEDIC displacement maps for all ${total} selected runs; raw phase preflight is not required."
  elif (( present > 0 )); then
    echo "ERROR: Incomplete MEDIC resume outputs (${complete}/${total} selected runs have nonempty displacement maps)." >&2
    echo "ERROR: Missing or empty expected map: $first_missing" >&2
    echo "ERROR: Refusing to switch partially processed MEDIC runs to TOPUP. Restart with START_FROM_MODULE=fieldmaps to complete or explicitly change distortion correction." >&2
    exit 2
  elif (( total > 0 )) && has_reusable_topup_outputs; then
    echo "Resuming with existing TOPUP fieldmaps (DISTORTION_CORRECTION_MODE=topup); no MEDIC displacement maps exist for the ${total} selected runs."
    DISTORTION_CORRECTION_MODE="topup"
    FUNC_NOFIELDMAP_MODE=0
  else
    echo "ERROR: Cannot resume ${START_FROM_MODULE}: no complete MEDIC displacement maps or reusable TOPUP outputs were found for ${total} selected runs." >&2
    echo "ERROR: Restore the required outputs or restart with START_FROM_MODULE=fieldmaps. Direct-B0 and no-fieldmap outputs require their explicit configuration modes." >&2
    exit 2
  fi
}

resolve_distortion_correction_mode() {
  local check_distortion_inputs=0
  local fieldmaps_idx stop_idx
  if should_run_stage "fieldmaps"; then
    check_distortion_inputs=1
    if [[ -n "${STOP_AFTER_MODULE:-}" ]]; then
      fieldmaps_idx="$(stage_index "fieldmaps")" || return 1
      stop_idx="$(stage_index "$STOP_AFTER_MODULE")" || return 1
      if (( stop_idx < fieldmaps_idx )); then
        check_distortion_inputs=0
      fi
    fi
  fi

  case "$DISTORTION_CORRECTION_MODE" in
    medic)
      if (( ! check_distortion_inputs )); then
        case "$START_FROM_MODULE" in
          coreg|headmotion) resolve_medic_resume_mode ;;
          *) echo "Skipping MEDIC raw-input preflight (START_FROM_MODULE=${START_FROM_MODULE})" ;;
        esac
      elif medic_phase_inputs_complete; then
        echo "MEDIC phase companions detected; using DISTORTION_CORRECTION_MODE=medic"
      elif has_topup_fieldmaps; then
        echo "WARNING: DISTORTION_CORRECTION_MODE=medic requires per-echo phase images, but complete phase companions were not found." >&2
        echo "WARNING: Falling back to DISTORTION_CORRECTION_MODE=topup because AP/PA fieldmaps are available." >&2
        DISTORTION_CORRECTION_MODE="topup"
      else
        echo "ERROR: DISTORTION_CORRECTION_MODE=medic requires per-echo phase images, and TOPUP fallback fieldmaps were not found." >&2
        echo "ERROR: Provide MEDIC phase companions or AP/PA fieldmaps. To run without distortion correction, set DISTORTION_CORRECTION_MODE=none; this can proceed but may produce suboptimal registration to anatomy." >&2
        exit 2
      fi
      ;;
    topup)
      if (( check_distortion_inputs )) && ! has_topup_fieldmaps; then
        echo "ERROR: DISTORTION_CORRECTION_MODE=topup requires AP/PA fieldmaps, but they were not found under func/unprocessed/${FUNC_DIRNAME}/field_maps." >&2
        echo "ERROR: Provide AP/PA fieldmaps, choose DISTORTION_CORRECTION_MODE=medic with phase companions, or explicitly set DISTORTION_CORRECTION_MODE=none. The none mode can proceed but may yield suboptimal registration to anatomy." >&2
        exit 2
      fi
      ;;
    none)
      FUNC_NOFIELDMAP_MODE=1
      echo "WARNING: DISTORTION_CORRECTION_MODE=none; no susceptibility distortion correction will be applied, which can degrade registration to anatomy." >&2
      ;;
  esac
}
: "${MGTR_INPUT_TAG:=}"
: "${MGTR_OUTPUT_TAG:=}"
: "${MGTR_ENABLE:=auto}" # auto|0|1; auto skips MGTR after aCompCor
: "${PIPELINE_SOURCE_FUNC_TAG:=}"

# Module-specific thread controls.
: "${THREADS_DEFAULT:=8}"
: "${THREADS_ANAT_HCP:=$THREADS_DEFAULT}"
: "${THREADS_FIELDMAPS:=$THREADS_DEFAULT}"
: "${THREADS_COREG:=$THREADS_DEFAULT}"
: "${THREADS_HEADMOTION:=$THREADS_DEFAULT}"
: "${THREADS_MEICA:=$THREADS_DEFAULT}"
: "${PIPELINE_QUIET_MODULE_OUTPUT:=1}" # 0|1
: "${PIPELINE_LOG_TAIL_LINES:=40}"

stage_index() {
  case "$1" in
    validate) echo 5 ;;
    anat_hcp) echo 10 ;;
    anat_charm) echo 20 ;;
    fieldmaps) echo 30 ;;
    coreg) echo 40 ;;
    headmotion) echo 50 ;;
    denoise|meica) echo 60 ;;
    mgtr) echo 70 ;;
    vol2surf) echo 80 ;;
    concat) echo 90 ;;
    nsi) echo 100 ;;
    pfm) echo 110 ;;
    pfm_update|infomap_update) echo 115 ;;
    tans) echo 120 ;;
    *)
      echo "ERROR: invalid module tag: $1" >&2
      return 1
      ;;
  esac
}

should_run_stage() {
  local stage="$1"
  local start_idx stage_idx
  start_idx="$(stage_index "$START_FROM_MODULE")" || return 1
  stage_idx="$(stage_index "$stage")" || return 1
  [ "$stage_idx" -ge "$start_idx" ]
}

is_all_task_selection() {
  local value="${1:-}"
  [[ "${value,,}" == "all" ]]
}

discover_task_dirs() {
  local root task_dir task
  declare -A seen=()
  for root in "$SubjectDir/func" "$SubjectDir/func/unprocessed"; do
    [[ -d "$root" ]] || continue
    while IFS= read -r task_dir; do
      task="$(basename "$task_dir")"
      [[ -n "$task" ]] || continue
      [[ -n "${seen[$task]:-}" ]] && continue
      compgen -G "$task_dir/session_*" > /dev/null || continue
      seen["$task"]=1
      printf '%s\n' "$task"
    done < <(find -L "$root" -mindepth 1 -maxdepth 1 -type d | sort -V)
  done
}

run_all_tasks() {
  local -a tasks=()
  local task rc=0
  mapfile -t tasks < <(discover_task_dirs)
  if [[ "${#tasks[@]}" -eq 0 ]]; then
    echo "ERROR: FUNC_DIRNAME=All but no task directories were found under $SubjectDir/func or $SubjectDir/func/unprocessed" >&2
    exit 2
  fi

  echo "Task selection: All"
  printf 'Discovered tasks (%d): %s\n' "${#tasks[@]}" "$(IFS=,; echo "${tasks[*]}")"

  local start_idx stop_idx
  local validate_idx anat_hcp_idx anat_charm_idx
  validate_idx="$(stage_index "validate")"
  anat_hcp_idx="$(stage_index "anat_hcp")"
  anat_charm_idx="$(stage_index "anat_charm")"

  start_idx="$(stage_index "$START_FROM_MODULE")" || exit 2
  if [[ -n "${STOP_AFTER_MODULE:-}" ]]; then
    stop_idx="$(stage_index "$STOP_AFTER_MODULE")" || exit 2
  else
    stop_idx=999999
  fi

  local needs_anat=0
  if (( start_idx <= anat_charm_idx && stop_idx >= anat_hcp_idx )); then
    needs_anat=1
  fi

  # Run shared (subject-level) anatomy once, then run per-task modules in child invocations.
  if (( needs_anat )); then
    local global_start="$START_FROM_MODULE"
    local global_stop="anat_charm"

    # Never run validate in the global/anat pass; validate is per-task.
    if (( start_idx <= validate_idx )); then
      global_start="anat_hcp"
    fi
    if (( start_idx > validate_idx && start_idx < anat_hcp_idx )); then
      global_start="anat_hcp"
    fi

    if (( stop_idx <= anat_charm_idx )); then
      global_stop="$STOP_AFTER_MODULE"
    fi

    echo
    echo "[task-loop] running anatomy once (start=${global_start} stop=${global_stop})"
    set +e
    env \
      PIPELINE_TASK_LOOP_CHILD=1 \
      PIPELINE_TASK_OVERRIDE="${tasks[0]}" \
      PIPELINE_TASK_XFMS_OVERRIDE="${tasks[0]}" \
      PIPELINE_SKIP_ANAT=0 \
      PIPELINE_START_FROM_MODULE_OVERRIDE="$global_start" \
      PIPELINE_STOP_AFTER_MODULE_OVERRIDE="$global_stop" \
      bash "$0" "$SubjectDir" "$ConfigFile"
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
      echo "[task-loop] FAILED anatomy pass exit=${rc}" >&2
      exit "$rc"
    fi

    if [[ -n "${STOP_AFTER_MODULE:-}" ]] && (( stop_idx <= anat_charm_idx )); then
      echo "[task-loop] stopping after ${STOP_AFTER_MODULE} (global stage); not running per-task modules"
      exit 0
    fi
  fi

  for task in "${tasks[@]}"; do
    echo
    echo "[task-loop] starting task=${task}"
    set +e
    env \
      PIPELINE_TASK_LOOP_CHILD=1 \
      PIPELINE_TASK_OVERRIDE="$task" \
      PIPELINE_TASK_XFMS_OVERRIDE="$task" \
      PIPELINE_SKIP_ANAT=1 \
      bash "$0" "$SubjectDir" "$ConfigFile"
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
      echo "[task-loop] FAILED task=${task} exit=${rc}" >&2
      exit "$rc"
    fi
    echo "[task-loop] complete task=${task}"
  done
  exit 0
}

if is_all_task_selection "$FUNC_DIRNAME" && [[ "$PIPELINE_TASK_LOOP_CHILD" != "1" ]]; then
  run_all_tasks
fi

# Apply optional post-config overrides.
if [[ -n "${PIPELINE_START_FROM_MODULE_OVERRIDE:-}" ]]; then
  START_FROM_MODULE="${PIPELINE_START_FROM_MODULE_OVERRIDE}"
fi
if [[ -n "${PIPELINE_STOP_AFTER_MODULE_OVERRIDE:-}" ]]; then
  STOP_AFTER_MODULE="${PIPELINE_STOP_AFTER_MODULE_OVERRIDE}"
fi
if [[ "${START_FROM_MODULE}" == "meica" ]]; then
  START_FROM_MODULE="denoise"
fi
if [[ "${STOP_AFTER_MODULE:-}" == "meica" ]]; then
  STOP_AFTER_MODULE="denoise"
fi

case "${TANS_ENABLE}" in
  0|1) ;;
  *) echo "ERROR: TANS_ENABLE must be 0 or 1" >&2; exit 2 ;;
esac
if [[ "${TANS_ENABLE}" != "1" ]] && \
   [[ "$START_FROM_MODULE" == "tans" || "$STOP_AFTER_MODULE" == "tans" ]]; then
  echo "ERROR: selecting the tans stage requires TANS_ENABLE=1" >&2
  exit 2
fi

PIPELINE_REQUESTED_DISTORTION_CORRECTION_MODE="$DISTORTION_CORRECTION_MODE"
resolve_distortion_correction_mode
export DISTORTION_CORRECTION_MODE FUNC_NOFIELDMAP_MODE AROMA_NSI_RECLASSIFY_ENABLE

if [ ! -f "$EnvironmentScript" ]; then
  echo "ERROR: missing environment setup script: $EnvironmentScript"
  exit 2
fi
source "$EnvironmentScript"

ensure_tkregister_compat() {
  if command -v tkregister >/dev/null 2>&1; then
    return 0
  fi
  local tkregister2_bin
  if ! tkregister2_bin="$(command -v tkregister2 2>/dev/null)"; then
    echo "ERROR: required command 'tkregister' is missing and fallback 'tkregister2' was not found." >&2
    return 1
  fi
  local shim_dir="${TMPDIR:-/tmp}/mefmri_compat_bin"
  mkdir -p "$shim_dir"
  cat > "${shim_dir}/tkregister" <<EOF
#!/usr/bin/env bash
exec "${tkregister2_bin}" "\$@"
EOF
  chmod +x "${shim_dir}/tkregister"
  export PATH="${shim_dir}:${PATH}"
  echo "INFO: Installed tkregister compatibility shim -> ${tkregister2_bin}"
}

ensure_tkregister_compat

resolve_freesurfer_license() {
  if [[ -n "${FS_LICENSE:-}" && -f "${FS_LICENSE}" ]]; then
    echo "${FS_LICENSE}"
    return 0
  fi
  if [[ -n "${FS_LICENSE_FILE:-}" && -f "${FS_LICENSE_FILE}" ]]; then
    echo "${FS_LICENSE_FILE}"
    return 0
  fi
  if [[ -n "${FREESURFER_HOME:-}" && -f "${FREESURFER_HOME}/license.txt" ]]; then
    echo "${FREESURFER_HOME}/license.txt"
    return 0
  fi
  return 1
}

# Anatomical HCP/FreeSurfer stages require a valid license file.
if should_run_stage "anat_hcp" && [[ "${PIPELINE_SKIP_ANAT}" != "1" ]]; then
  if ! FS_LICENSE_EFFECTIVE="$(resolve_freesurfer_license)"; then
    echo "ERROR: FreeSurfer license not found." >&2
    echo "Set FS_LICENSE or FS_LICENSE_FILE to a readable license.txt path." >&2
    exit 2
  fi
  export FS_LICENSE="${FS_LICENSE_EFFECTIVE}"
fi

case "${AtlasSpace}" in
  T1w|t1w|Tlw|tlw) AtlasSpace="T1w" ;;
  MNINonlinear|mni|MNI|mninonlinear) AtlasSpace="MNINonlinear" ;;
  *)
    echo "ERROR: invalid AtlasSpace='$AtlasSpace' (supported: T1w or MNINonlinear)"
    exit 2
    ;;
esac

case "${PIPELINE_QUIET_MODULE_OUTPUT}" in
  0|1) ;;
  *)
    echo "ERROR: invalid PIPELINE_QUIET_MODULE_OUTPUT='${PIPELINE_QUIET_MODULE_OUTPUT}' (expected 0 or 1)"
    exit 2
    ;;
esac

case "${COREG_SDC_QC_ENABLE}" in
  0|1) ;;
  *)
    echo "ERROR: invalid COREG_SDC_QC_ENABLE='${COREG_SDC_QC_ENABLE}' (expected 0 or 1)"
    exit 2
    ;;
esac

case "${QC_VALIDATION_ENABLE}" in
  0|1) ;;
  *)
    echo "ERROR: invalid QC_VALIDATION_ENABLE='${QC_VALIDATION_ENABLE}' (expected 0 or 1)"
    exit 2
    ;;
esac

case "${COREG_SDC_QC_STRICT}" in
  0|1) ;;
  *)
    echo "ERROR: invalid COREG_SDC_QC_STRICT='${COREG_SDC_QC_STRICT}' (expected 0 or 1)"
    exit 2
    ;;
esac

case "${MOTION_GRAYPLOTS_QC_ENABLE}" in
  0|1) ;;
  *)
    echo "ERROR: invalid MOTION_GRAYPLOTS_QC_ENABLE='${MOTION_GRAYPLOTS_QC_ENABLE}' (expected 0 or 1)"
    exit 2
    ;;
esac

case "${MOTION_GRAYPLOTS_QC_STRICT}" in
  0|1) ;;
  *)
    echo "ERROR: invalid MOTION_GRAYPLOTS_QC_STRICT='${MOTION_GRAYPLOTS_QC_STRICT}' (expected 0 or 1)"
    exit 2
    ;;
esac

if ! [[ "${MOTION_GRAYPLOTS_QC_FD_THRESHOLD}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: invalid MOTION_GRAYPLOTS_QC_FD_THRESHOLD='${MOTION_GRAYPLOTS_QC_FD_THRESHOLD}' (expected non-negative numeric mm threshold)"
  exit 2
fi

case "${PROCESSING_MODE}" in
  auto|multi_echo|single_echo) ;;
  *)
    echo "ERROR: invalid PROCESSING_MODE='${PROCESSING_MODE}' (expected auto, multi_echo, or single_echo)"
    exit 2
    ;;
esac

case "${SINGLE_ECHO_DENOISE_METHOD}" in
  aroma|acompcor) ;;
  *)
    echo "ERROR: invalid SINGLE_ECHO_DENOISE_METHOD='${SINGLE_ECHO_DENOISE_METHOD}' (expected aroma or acompcor)"
    exit 2
    ;;
esac

case "${MULTI_ECHO_DENOISE_METHOD}" in
  meica|aroma|acompcor) ;;
  *)
    echo "ERROR: invalid MULTI_ECHO_DENOISE_METHOD='${MULTI_ECHO_DENOISE_METHOD}' (expected meica, aroma, or acompcor)"
    exit 2
    ;;
esac

if ! [[ "${SINGLE_ECHO_ECHO_INDEX}" =~ ^[0-9]+$ ]] || [[ "${SINGLE_ECHO_ECHO_INDEX}" -lt 1 ]]; then
  echo "ERROR: SINGLE_ECHO_ECHO_INDEX must be an integer >= 1 (got '${SINGLE_ECHO_ECHO_INDEX}')"
  exit 2
fi

detect_run_echo_count() {
  local run_dir="$1"
  local te_file
  for te_file in "$run_dir/TE.txt" "$run_dir/te.txt"; do
    if [[ -f "$te_file" ]]; then
      awk 'NF{print NF; exit}' "$te_file"
      return 0
    fi
  done

  local count=0
  shopt -s nullglob
  local matches=( "$run_dir"/"${FUNC_FILE_PREFIX}"*_E*.nii.gz )
  shopt -u nullglob
  local f base
  count=0
  for f in "${matches[@]}"; do
    base="$(basename "$f")"
    if [[ "$base" =~ ^${FUNC_FILE_PREFIX}.*_E[0-9]+(_acpc)?\.nii\.gz$ ]]; then
      count=$((count + 1))
    fi
  done
  if [[ "$count" -gt 0 ]]; then
    echo "$count"
    return 0
  fi

  echo ""
}

detect_dataset_min_echoes() {
  local -a roots=(
    "$SubjectDir/func/$FUNC_DIRNAME"
    "$SubjectDir/func/unprocessed/$FUNC_DIRNAME"
  )
  local root run_dir count min_count=""
  for root in "${roots[@]}"; do
    [[ -d "$root" ]] || continue
    while IFS= read -r run_dir; do
      count="$(detect_run_echo_count "$run_dir")"
      [[ -n "$count" ]] || continue
      if [[ -z "$min_count" || "$count" -lt "$min_count" ]]; then
        min_count="$count"
      fi
    done < <(find -L "$root" -mindepth 2 -maxdepth 2 -type d -name 'run_*' | sort -V)
  done
  if [[ -z "$min_count" ]]; then
    echo "0"
  else
    echo "$min_count"
  fi
}

PIPELINE_MIN_ECHOES="$(detect_dataset_min_echoes)"
PIPELINE_DENOISE_FALLBACK_REASON=""
PIPELINE_EFFECTIVE_DENOISE_MODE="multi_echo"
if [[ "${PROCESSING_MODE}" == "single_echo" ]]; then
  PIPELINE_EFFECTIVE_DENOISE_MODE="single_echo"
  PIPELINE_DENOISE_FALLBACK_REASON="explicit_single_echo"
elif [[ "${PIPELINE_MIN_ECHOES}" -lt 3 ]]; then
  PIPELINE_EFFECTIVE_DENOISE_MODE="single_echo"
  PIPELINE_DENOISE_FALLBACK_REASON="echo_count_lt_3"
fi

PIPELINE_SOURCE_FUNC_TAG="OCME"
PIPELINE_EFFECTIVE_DENOISE_METHOD="$MULTI_ECHO_DENOISE_METHOD"
PIPELINE_DENOISE_OUTPUT_TAG="OCME+MEICA"
if [[ "${PIPELINE_EFFECTIVE_DENOISE_MODE}" == "single_echo" ]]; then
  PIPELINE_SOURCE_FUNC_TAG="E${SINGLE_ECHO_ECHO_INDEX}"
  PIPELINE_EFFECTIVE_DENOISE_METHOD="$SINGLE_ECHO_DENOISE_METHOD"
  if [[ "${SINGLE_ECHO_DENOISE_METHOD}" == "aroma" ]]; then
    PIPELINE_DENOISE_OUTPUT_TAG="${PIPELINE_SOURCE_FUNC_TAG}+AROMA"
  else
    PIPELINE_DENOISE_OUTPUT_TAG="${PIPELINE_SOURCE_FUNC_TAG}+aCompCor"
  fi
else
  case "${MULTI_ECHO_DENOISE_METHOD}" in
    meica)
      PIPELINE_DENOISE_OUTPUT_TAG="OCME+MEICA"
      ;;
    aroma)
      PIPELINE_DENOISE_OUTPUT_TAG="OCME+AROMA"
      ;;
    acompcor)
      PIPELINE_DENOISE_OUTPUT_TAG="OCME+aCompCor"
      ;;
  esac
fi
PIPELINE_MGTR_OUTPUT_TAG_DEFAULT="${PIPELINE_DENOISE_OUTPUT_TAG}+MGTR"
PIPELINE_MGTR_TAGS_CONFIGURED=0
if [[ -n "${MGTR_INPUT_TAG}" || -n "${MGTR_OUTPUT_TAG}" ]]; then
  PIPELINE_MGTR_TAGS_CONFIGURED=1
fi
case "${MGTR_ENABLE}" in
  auto)
    if [[ "${PIPELINE_EFFECTIVE_DENOISE_METHOD}" == "acompcor" && "${PIPELINE_MGTR_TAGS_CONFIGURED}" == "0" ]]; then
      PIPELINE_EFFECTIVE_MGTR_ENABLE=0
    else
      PIPELINE_EFFECTIVE_MGTR_ENABLE=1
    fi
    ;;
  0|1)
    PIPELINE_EFFECTIVE_MGTR_ENABLE="$MGTR_ENABLE"
    ;;
  *)
    echo "ERROR: invalid MGTR_ENABLE='${MGTR_ENABLE}' (expected auto, 0, or 1)"
    exit 2
    ;;
esac
: "${MGTR_INPUT_TAG:=$PIPELINE_DENOISE_OUTPUT_TAG}"
: "${MGTR_OUTPUT_TAG:=$PIPELINE_MGTR_OUTPUT_TAG_DEFAULT}"
if [[ -z "${VOL2SURF_INPUTS}" ]]; then
  if [[ "${PIPELINE_EFFECTIVE_MGTR_ENABLE}" == "1" ]]; then
    VOL2SURF_INPUTS="${PIPELINE_SOURCE_FUNC_TAG},${PIPELINE_DENOISE_OUTPUT_TAG},${MGTR_OUTPUT_TAG}"
  else
    VOL2SURF_INPUTS="${PIPELINE_SOURCE_FUNC_TAG},${PIPELINE_DENOISE_OUTPUT_TAG}"
  fi
fi
if [[ "${PIPELINE_EFFECTIVE_MGTR_ENABLE}" == "1" ]]; then
  : "${CONCAT_INPUT_TAG:=$MGTR_OUTPUT_TAG}"
else
  : "${CONCAT_INPUT_TAG:=$PIPELINE_DENOISE_OUTPUT_TAG}"
fi
: "${NSI_INPUT_TAG:=$CONCAT_INPUT_TAG}"
: "${PFM_INPUT_TAG:=$CONCAT_INPUT_TAG}"

export PIPELINE_MIN_ECHOES
export PIPELINE_EFFECTIVE_DENOISE_MODE
export PIPELINE_DENOISE_FALLBACK_REASON
export PIPELINE_DENOISE_OUTPUT_TAG
export PIPELINE_SOURCE_FUNC_TAG
export PIPELINE_EFFECTIVE_DENOISE_METHOD
export PIPELINE_EFFECTIVE_MGTR_ENABLE
export MGTR_INPUT_TAG
export MGTR_OUTPUT_TAG

echo
echo "ME-fMRI Pipeline"
echo "MEDIR: ${MEDIR}"
echo "SubjectDir: ${SubjectDir}"
echo "Subject: ${Subject}"
echo "StudyFolder: ${StudyFolder}"
echo "START_SESSION: ${START_SESSION}"
echo "START_FROM_MODULE: ${START_FROM_MODULE}"
echo "STOP_AFTER_MODULE: ${STOP_AFTER_MODULE:-<none>}"
echo "AtlasSpace: ${AtlasSpace}"
echo "Task loop child: ${PIPELINE_TASK_LOOP_CHILD}"
echo "Skip anatomy: ${PIPELINE_SKIP_ANAT}"
echo "Functional naming: func/${FUNC_DIRNAME}, prefix ${FUNC_FILE_PREFIX}_*"
echo "Functional xfm namespace: func/xfms/${FUNC_XFMS_DIRNAME}"
echo "Distortion correction mode: ${DISTORTION_CORRECTION_MODE}"
echo "FUNC_NOFIELDMAP_MODE: ${FUNC_NOFIELDMAP_MODE}"
echo "MGTR: requested=${MGTR_ENABLE} effective=${PIPELINE_EFFECTIVE_MGTR_ENABLE}"
echo "Processing mode: requested=${PROCESSING_MODE} effective=${PIPELINE_EFFECTIVE_DENOISE_MODE} min_echoes=${PIPELINE_MIN_ECHOES}"
if [[ "${PIPELINE_EFFECTIVE_DENOISE_MODE}" == "single_echo" ]]; then
  echo "Single-echo denoise method: ${SINGLE_ECHO_DENOISE_METHOD} (source echo E${SINGLE_ECHO_ECHO_INDEX})"
  if [[ -n "${PIPELINE_DENOISE_FALLBACK_REASON}" ]]; then
    echo "Single-echo selection reason: ${PIPELINE_DENOISE_FALLBACK_REASON}"
  fi
else
  echo "Multi-echo denoise method: ${MULTI_ECHO_DENOISE_METHOD}"
fi
echo "Threads (anat_hcp,fieldmaps,coreg,headmotion,denoise): ${THREADS_ANAT_HCP},${THREADS_FIELDMAPS},${THREADS_COREG},${THREADS_HEADMOTION},${THREADS_MEICA}"
echo "MEICA defaults: tedana_env=${TEDANA_ENV:-unset}, compat=${TEDANA_COMPAT_MODE:-unset}, pca=${MEPCA}"
echo "Masking defaults: CHARM_BRAIN_MASK_MODE=${CHARM_BRAIN_MASK_MODE:-unset}, VOL2SURF_USE_CORTICAL_RIBBON_MASK=${VOL2SURF_USE_CORTICAL_RIBBON_MASK:-unset}, VOL2SURF_USE_GOOD_VOXELS_MASK=${VOL2SURF_USE_GOOD_VOXELS_MASK:-${VOL2SURF_USE_GOOD_VOXELS_MASK_FINAL:-unset}}"
echo "Reclass defaults: mode=${MEICA_CLASSIFIER_MODE:-unset}, nsi_kill=${MEICA_NSI_KILL_MODE:-unset}, reports_disabled=${MEICA_RECLASS_NO_REPORTS:-unset}"
echo "Module output mode: quiet=${PIPELINE_QUIET_MODULE_OUTPUT}"

FUNC_QC_DIR="$SubjectDir/func/$FUNC_DIRNAME/$QC_OUTPUT_DIRNAME"
FUNC_QA_DIR="$FUNC_QC_DIR"
MODULE_LOG_DIR="$FUNC_QC_DIR/ModuleLogs"
mkdir -p "$MODULE_LOG_DIR"

run_module() {
  local module="$1"
  shift

  if [[ "${PIPELINE_QUIET_MODULE_OUTPUT}" == "0" ]]; then
    "$@"
    return $?
  fi

  # Some denoising methods share the same pipeline scheduling stage.
  # Allow callers to override the user-facing label and logfile suffix without
  # changing scheduling semantics.
  local display_tag log_tag
  display_tag="${MODULE_DISPLAY_TAG:-$module}"
  log_tag="${MODULE_LOG_TAG:-$module}"

  local ts logfile rc
  ts="$(date +%Y%m%d_%H%M%S)"
  logfile="$MODULE_LOG_DIR/${ts}_${log_tag}.log"
  echo "[$display_tag] log: $logfile"

  # Always capture complete output to log; only stream high-level markers to terminal.
  set +e
  "$@" 2>&1 | awk -v logfile_path="$logfile" '
    {
      print >> logfile_path
      fflush(logfile_path)
      if ($0 ~ /^\[(concat|crawling_seed_fc|nsi|pfm)\] complete$/) {
        next
      }
      if ($0 ~ /^\[coreg\]/ || $0 ~ /^\[headmotion\]/ || $0 ~ /^\[medic\]/ || $0 ~ /^\[concat\]/ || $0 ~ /^\[nsi\]/ || $0 ~ /^\[pfm\]/ ||
          $0 ~ /^\[crawling_seed_fc\]/ ||
          $0 ~ /^\[qc-coreg-sdc\]/ ||
          $0 ~ /^\[qc-motion-grayplots\]/ ||
          $0 ~ /^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\] MEICA: start subject=/ ||
          $0 ~ /^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\] MEICA: start session_/ ||
          $0 ~ /^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\] MEICA: done session_/ ||
          $0 ~ /^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\] MEICA: all done subject=/) {
        print
      }
    }
  '
  rc=${PIPESTATUS[0]}
  set -e

  if [[ "$rc" -ne 0 ]]; then
    echo "[$display_tag] FAILED (exit $rc). Log: $logfile" >&2
    echo "[$display_tag] error summary:" >&2
    if ! tail -n 400 "$logfile" | awk '
      /Traceback \(most recent call last\):/ ||
      /^ERROR:/ ||
      /FileNotFoundError:/ ||
      /Exception:/ ||
      /FAILED/ ||
      /can'\''t open file/ { print; found=1 }
      END { exit(found?0:1) }
    ' >&2; then
      echo "[$display_tag] (no explicit traceback/error markers found; showing last ${PIPELINE_LOG_TAIL_LINES} lines)" >&2
      tail -n "$PIPELINE_LOG_TAIL_LINES" "$logfile" >&2 || true
    fi
    return "$rc"
  fi

  echo "[$display_tag] complete"
}

print_section() {
  local title="$1"
  echo
  echo "$title"
}

run_motion_grayplots_qc() {
  if [[ "${MOTION_GRAYPLOTS_QC_ENABLE}" != "1" ]]; then
    echo "Skipping Motion & GrayPlots QC (MOTION_GRAYPLOTS_QC_ENABLE=${MOTION_GRAYPLOTS_QC_ENABLE})"
    return 0
  fi
  [ -f "$FUNC_QC_MOTION_GRAYPLOTS_MODULE" ] || { echo "ERROR: missing module: $FUNC_QC_MOTION_GRAYPLOTS_MODULE"; exit 2; }
  run_module "qc_motion_grayplots" bash "$FUNC_QC_MOTION_GRAYPLOTS_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
}

if [[ "${RUN_CONFIG_SNAPSHOT}" == "1" ]]; then
  RUN_META_DIR="$FUNC_QC_DIR/RunMetadata"
  mkdir -p "$RUN_META_DIR"
  RUN_META_FILE="$RUN_META_DIR/pipeline_run_$(date +%Y%m%d_%H%M%S).txt"
  {
    echo "date=$(date --iso-8601=seconds)"
    echo "subject_dir=$SubjectDir"
    echo "config_file=$ConfigFile"
    echo "start_from_module=$START_FROM_MODULE"
    echo "stop_after_module=${STOP_AFTER_MODULE:-}"
    echo "atlas_space=$AtlasSpace"
    echo "task_loop_child=$PIPELINE_TASK_LOOP_CHILD"
    echo "func_dirname=$FUNC_DIRNAME"
    echo "func_file_prefix=$FUNC_FILE_PREFIX"
    echo "func_xfms_dirname=$FUNC_XFMS_DIRNAME"
    echo "distortion_correction_mode_requested=$PIPELINE_REQUESTED_DISTORTION_CORRECTION_MODE"
    echo "distortion_correction_mode=$DISTORTION_CORRECTION_MODE"
    echo "qc_output_dirname=$QC_OUTPUT_DIRNAME"
    echo "qc_validation_enable=$QC_VALIDATION_ENABLE"
    echo "coreg_sdc_qc_enable=$COREG_SDC_QC_ENABLE"
    echo "coreg_sdc_qc_strict=$COREG_SDC_QC_STRICT"
    echo "motion_grayplots_qc_enable=$MOTION_GRAYPLOTS_QC_ENABLE"
    echo "motion_grayplots_qc_strict=$MOTION_GRAYPLOTS_QC_STRICT"
    echo "motion_grayplots_qc_fd_threshold=$MOTION_GRAYPLOTS_QC_FD_THRESHOLD"
    echo "sbref_echo_combination=${SBREF_ECHO_COMBINATION:-}"
    echo "sbref_max_te_ms=${SBREF_MAX_TE_MS:-}"
    echo "sbref_t2smap_mask_thr=${SBREF_T2SMAP_MASK_THR:-}"
    echo "sbref_t2smap_threads=${SBREF_T2SMAP_THREADS:-}"
    echo "sbref_keep_intermediates=${SBREF_KEEP_INTERMEDIATES:-}"
    echo "processing_mode_requested=${PROCESSING_MODE}"
    echo "processing_mode_effective=${PIPELINE_EFFECTIVE_DENOISE_MODE}"
    echo "processing_mode_reason=${PIPELINE_DENOISE_FALLBACK_REASON}"
    echo "min_echoes=${PIPELINE_MIN_ECHOES}"
    echo "single_echo_method=${SINGLE_ECHO_DENOISE_METHOD}"
    echo "multi_echo_method=${MULTI_ECHO_DENOISE_METHOD}"
    echo "single_echo_echo_index=${SINGLE_ECHO_ECHO_INDEX}"
    echo "aroma_nsi_reclassify_enable=${AROMA_NSI_RECLASSIFY_ENABLE}"
    echo "aroma_nsi_threshold=${AROMA_NSI_THRESHOLD}"
    echo "pipeline_source_func_tag=${PIPELINE_SOURCE_FUNC_TAG}"
    echo "pipeline_denoise_output_tag=${PIPELINE_DENOISE_OUTPUT_TAG}"
    echo "pipeline_effective_denoise_method=${PIPELINE_EFFECTIVE_DENOISE_METHOD}"
    echo "mgtr_enable=${MGTR_ENABLE}"
    echo "pipeline_effective_mgtr_enable=${PIPELINE_EFFECTIVE_MGTR_ENABLE}"
    echo "mgtr_input_tag=${MGTR_INPUT_TAG}"
    echo "mgtr_output_tag=${MGTR_OUTPUT_TAG}"
    echo "concat_input_tag=${CONCAT_INPUT_TAG:-}"
    echo "nsi_input_tag=${NSI_INPUT_TAG:-}"
    echo "pfm_input_tag=${PFM_INPUT_TAG:-}"
    echo "charm_brain_mask_mode=${CHARM_BRAIN_MASK_MODE:-}"
    echo "vol2surf_use_cortical_ribbon_mask=${VOL2SURF_USE_CORTICAL_RIBBON_MASK:-}"
    echo "vol2surf_use_good_voxels_mask=${VOL2SURF_USE_GOOD_VOXELS_MASK:-${VOL2SURF_USE_GOOD_VOXELS_MASK_FINAL:-}}"
    echo "vol2surf_good_voxels_factor=${VOL2SURF_GOOD_VOXELS_FACTOR:-}"
    echo "vol2surf_good_voxels_sigma_mm=${VOL2SURF_GOOD_VOXELS_SIGMA_MM:-}"
    echo "vol2surf_good_voxels_keep_intermediates=${VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES:-}"
    echo "tedana_env=${TEDANA_ENV:-}"
    echo "tedana_compat_mode=${TEDANA_COMPAT_MODE:-}"
    echo "warpkit_env=${WARPKIT_ENV:-}"
    echo "warpkit_activate_mode=${WARPKIT_ACTIVATE_MODE:-}"
    echo "mepca=$MEPCA"
    echo "meica_classifier_mode=${MEICA_CLASSIFIER_MODE:-}"
    echo "meica_nsi_kill_mode=${MEICA_NSI_KILL_MODE:-}"
    echo "meica_reclass_no_reports=${MEICA_RECLASS_NO_REPORTS:-}"
    echo "concat_enable=${CONCAT_ENABLE}"
    echo "nsi_enable=${NSI_ENABLE}"
    echo "pfm_enable=${PFM_ENABLE}"
    echo "pfm_infomap_update_enable=${PFM_INFOMAP_UPDATE_ENABLE}"
    echo "crawling_seed_fc_enable=${CRAWLING_SEED_FC_ENABLE}"
  } > "$RUN_META_FILE"
  echo "Run metadata snapshot: $RUN_META_FILE"
fi

print_section "Running input validation"
if should_run_stage "validate"; then
  [ -f "$VALIDATE_MODULE" ] || { echo "ERROR: missing module: $VALIDATE_MODULE"; exit 2; }
  run_module "validate" bash "$VALIDATE_MODULE" "$SubjectDir" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX" "$START_SESSION"
  if [[ "$STOP_AFTER_MODULE" == "validate" ]]; then
    echo "Stopping after validate (STOP_AFTER_MODULE=validate)"
    exit 0
  fi
else
  echo "Skipping validate (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Running anatomical modules"
if should_run_stage "anat_hcp"; then
  if [[ "${PIPELINE_SKIP_ANAT}" == "1" ]]; then
    echo "Skipping anat_hcp (PIPELINE_SKIP_ANAT=1)"
  else
    [ -f "$ANAT_HCP_MODULE" ] || { echo "ERROR: missing module: $ANAT_HCP_MODULE"; exit 2; }
    run_module "anat_hcp" bash "$ANAT_HCP_MODULE" "$StudyFolder" "$Subject" "$THREADS_ANAT_HCP"
  fi
  if [[ "$STOP_AFTER_MODULE" == "anat_hcp" ]]; then
    echo "Stopping after anat_hcp (STOP_AFTER_MODULE=anat_hcp)"
    exit 0
  fi
else
  echo "Skipping anat_hcp (START_FROM_MODULE=${START_FROM_MODULE})"
fi
if should_run_stage "anat_charm"; then
  if [[ "${PIPELINE_SKIP_ANAT}" == "1" ]]; then
    echo "Skipping anat_charm (PIPELINE_SKIP_ANAT=1)"
  else
    [ -f "$ANAT_CHARM_MODULE" ] || { echo "ERROR: missing module: $ANAT_CHARM_MODULE"; exit 2; }
    if [ -n "$CHARM_BIN" ]; then
      run_module "anat_charm" bash "$ANAT_CHARM_MODULE" "$StudyFolder" "$Subject" "$CHARM_BIN"
    else
      run_module "anat_charm" bash "$ANAT_CHARM_MODULE" "$StudyFolder" "$Subject"
    fi
  fi
  if [[ "$STOP_AFTER_MODULE" == "anat_charm" ]]; then
    echo "Stopping after anat_charm (STOP_AFTER_MODULE=anat_charm)"
    exit 0
  fi
else
  echo "Skipping anat_charm (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Processing distortion correction inputs"
if should_run_stage "fieldmaps"; then
  case "$DISTORTION_CORRECTION_MODE" in
    topup)
      [ -f "$FUNC_FIELDMAPS_MODULE" ] || { echo "ERROR: missing module: $FUNC_FIELDMAPS_MODULE"; exit 2; }
      [ -x "$FUNC_FIELDMAPS_MODULE" ] || chmod +x "$FUNC_FIELDMAPS_MODULE"
      run_module "fieldmaps" env FIELDMAP_STRATEGY=topup bash "$FUNC_FIELDMAPS_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$THREADS_FIELDMAPS" "$START_SESSION" "$FUNC_DIRNAME"
      ;;
    direct_b0|fieldmap|gre_fieldmap)
      [ -f "$FUNC_FIELDMAPS_MODULE" ] || { echo "ERROR: missing module: $FUNC_FIELDMAPS_MODULE"; exit 2; }
      [ -x "$FUNC_FIELDMAPS_MODULE" ] || chmod +x "$FUNC_FIELDMAPS_MODULE"
      run_module "fieldmaps" env FIELDMAP_STRATEGY=direct_b0 bash "$FUNC_FIELDMAPS_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$THREADS_FIELDMAPS" "$START_SESSION" "$FUNC_DIRNAME"
      ;;
    none)
      echo "Skipping fieldmap estimation (DISTORTION_CORRECTION_MODE=none; FUNC_NOFIELDMAP_MODE=1)"
      ;;
    medic)
      [ -f "$FUNC_MEDIC_MODULE" ] || { echo "ERROR: missing module: $FUNC_MEDIC_MODULE"; exit 2; }
      [ -x "$FUNC_MEDIC_MODULE" ] || chmod +x "$FUNC_MEDIC_MODULE"
      run_module "medic" env \
        WARPKIT_ENV="${WARPKIT_ENV:-warpkit_env}" \
        WARPKIT_ACTIVATE_MODE="${WARPKIT_ACTIVATE_MODE:-conda_activate}" \
        WARPKIT_MEDIC_BIN="${WARPKIT_MEDIC_BIN:-wk-medic}" \
        WARPKIT_COMPUTE_JACOBIAN_BIN="${WARPKIT_COMPUTE_JACOBIAN_BIN:-wk-compute-jacobian}" \
        MEDIC_OVERWRITE="${MEDIC_OVERWRITE:-0}" \
        MEDIC_NOISE_FRAMES="${MEDIC_NOISE_FRAMES:-0}" \
        MEDIC_DEBUG="${MEDIC_DEBUG:-0}" \
        MEDIC_WRAP_LIMIT="${MEDIC_WRAP_LIMIT:-0}" \
        MEDIC_COMPUTE_JACOBIAN="${MEDIC_COMPUTE_JACOBIAN:-1}" \
        FUNC_FILE_PREFIX="${FUNC_FILE_PREFIX:-Rest}" \
        bash "$FUNC_MEDIC_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$THREADS_FIELDMAPS" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
      ;;
  esac
  if [[ "$STOP_AFTER_MODULE" == "fieldmaps" ]]; then
    echo "Stopping after distortion correction inputs (STOP_AFTER_MODULE=fieldmaps)"
    exit 0
  fi
else
  echo "Skipping distortion correction inputs (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Coregistering SBrefs to anatomical image"
CoregModule="$FUNC_COREG_MODULE"
CoregTag="coreg"
if [[ "$DISTORTION_CORRECTION_MODE" == "medic" ]]; then
  CoregModule="$FUNC_MEDIC_COREG_MODULE"
  CoregTag="coreg_medic"
fi
[ -f "$CoregModule" ] || { echo "ERROR: missing module: $CoregModule"; exit 2; }
if should_run_stage "coreg"; then
  run_module "$CoregTag" bash "$CoregModule" "$MEDIR" "$Subject" "$StudyFolder" "$AtlasTemplate" "$DOF" "$THREADS_COREG" "$START_SESSION" "$AtlasSpace" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
  if [[ "${COREG_SDC_QC_ENABLE}" == "1" ]]; then
    [ -f "$FUNC_QC_COREG_SDC_MODULE" ] || { echo "ERROR: missing module: $FUNC_QC_COREG_SDC_MODULE"; exit 2; }
    run_module "qc_coreg_sdc" bash "$FUNC_QC_COREG_SDC_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX" "$FUNC_XFMS_DIRNAME" "$DISTORTION_CORRECTION_MODE"
  else
    echo "Skipping coreg/SDC QC report (COREG_SDC_QC_ENABLE=${COREG_SDC_QC_ENABLE})"
  fi
  if [[ "$STOP_AFTER_MODULE" == "coreg" ]]; then
    echo "Stopping after coreg (STOP_AFTER_MODULE=coreg)"
    exit 0
  fi
else
  echo "Skipping coreg (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Applying slice-time/headmotion/distortion corrections"
HeadmotionModule="$FUNC_HEADMOTION_MODULE"
HeadmotionTag="headmotion"
if [[ "$DISTORTION_CORRECTION_MODE" == "medic" ]]; then
  HeadmotionModule="$FUNC_MEDIC_HEADMOTION_MODULE"
  HeadmotionTag="headmotion_medic"
fi
[ -f "$HeadmotionModule" ] || { echo "ERROR: missing module: $HeadmotionModule"; exit 2; }
if should_run_stage "headmotion"; then
  run_module "$HeadmotionTag" bash "$HeadmotionModule" "$MEDIR" "$Subject" "$StudyFolder" "$AtlasTemplate" "$DOF" "$THREADS_HEADMOTION" "$START_SESSION" "$AtlasSpace" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
  if [[ "$STOP_AFTER_MODULE" == "headmotion" ]]; then
    echo "Stopping after headmotion (STOP_AFTER_MODULE=headmotion)"
    exit 0
  fi
else
  echo "Skipping headmotion (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Running denoising"
if should_run_stage "denoise"; then
  if [[ "${PIPELINE_EFFECTIVE_DENOISE_MODE}" == "single_echo" ]]; then
    if [[ "${SINGLE_ECHO_DENOISE_METHOD}" == "aroma" ]]; then
      [ -f "$FUNC_AROMA_MODULE" ] || { echo "ERROR: missing module: $FUNC_AROMA_MODULE"; exit 2; }
      MODULE_DISPLAY_TAG="${SINGLE_ECHO_DENOISE_METHOD}" MODULE_LOG_TAG="${SINGLE_ECHO_DENOISE_METHOD}" \
        run_module "denoise" bash "$FUNC_AROMA_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$START_SESSION" "$MEDIR"
    else
      [ -f "$FUNC_SINGLEECHO_MODULE" ] || { echo "ERROR: missing module: $FUNC_SINGLEECHO_MODULE"; exit 2; }
      MODULE_DISPLAY_TAG="${SINGLE_ECHO_DENOISE_METHOD}" MODULE_LOG_TAG="${SINGLE_ECHO_DENOISE_METHOD}" \
        run_module "denoise" bash "$FUNC_SINGLEECHO_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$START_SESSION" "$MEDIR"
    fi
  else
    [ -f "$FUNC_MEICA_MODULE" ] || { echo "ERROR: missing module: $FUNC_MEICA_MODULE"; exit 2; }
    if [[ "${MULTI_ECHO_DENOISE_METHOD}" == "meica" ]]; then
      run_module "denoise" bash "$FUNC_MEICA_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$MEPCA" "$MaxIterations" "$MaxRestarts" "$START_SESSION" "$MEDIR"
    else
      run_module "meica_optcom" env MEICA_RECLASSIFY_ENABLE=0 bash "$FUNC_MEICA_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$MEPCA" "$MaxIterations" "$MaxRestarts" "$START_SESSION" "$MEDIR"
      if [[ "${MULTI_ECHO_DENOISE_METHOD}" == "aroma" ]]; then
        [ -f "$FUNC_AROMA_MODULE" ] || { echo "ERROR: missing module: $FUNC_AROMA_MODULE"; exit 2; }
        MODULE_DISPLAY_TAG="${MULTI_ECHO_DENOISE_METHOD}" MODULE_LOG_TAG="${MULTI_ECHO_DENOISE_METHOD}" \
          run_module "denoise" bash "$FUNC_AROMA_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$START_SESSION" "$MEDIR"
      else
        [ -f "$FUNC_SINGLEECHO_MODULE" ] || { echo "ERROR: missing module: $FUNC_SINGLEECHO_MODULE"; exit 2; }
        MODULE_DISPLAY_TAG="${MULTI_ECHO_DENOISE_METHOD}" MODULE_LOG_TAG="${MULTI_ECHO_DENOISE_METHOD}" \
          run_module "denoise" bash "$FUNC_SINGLEECHO_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$START_SESSION" "$MEDIR"
      fi
    fi
  fi
  if [[ "$STOP_AFTER_MODULE" == "denoise" ]]; then
    echo "Skipping Motion & GrayPlots QC (requires vol2surf dtseries outputs; run through vol2surf or later)."
    echo "Stopping after denoise (STOP_AFTER_MODULE=denoise)"
    exit 0
  fi
else
  echo "Skipping denoise (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Running MGTR"
[ -f "$FUNC_MGTR_MODULE" ] || { echo "ERROR: missing module: $FUNC_MGTR_MODULE"; exit 2; }
if [[ "${PIPELINE_EFFECTIVE_MGTR_ENABLE}" != "1" ]]; then
  echo "Skipping mgtr (MGTR_ENABLE=${MGTR_ENABLE}, effective denoise method=${PIPELINE_EFFECTIVE_DENOISE_METHOD})"
elif should_run_stage "mgtr"; then
  run_module "mgtr" bash "$FUNC_MGTR_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION"
else
  echo "Skipping mgtr (START_FROM_MODULE=${START_FROM_MODULE})"
fi
if [[ "$STOP_AFTER_MODULE" == "mgtr" ]]; then
  echo "Skipping Motion & GrayPlots QC (requires vol2surf dtseries outputs; run through vol2surf or later)."
  echo "Stopping after mgtr (STOP_AFTER_MODULE=mgtr)"
  exit 0
fi

print_section "Mapping denoised data to surface"
[ -f "$FUNC_VOL2SURF_MODULE" ] || { echo "ERROR: missing module: $FUNC_VOL2SURF_MODULE"; exit 2; }
if [ -z "${VOL2SURF_INPUTS:-}" ]; then
  echo "ERROR: VOL2SURF_INPUTS is empty. Set a comma-separated list (e.g., OCME,OCME+MEICA,OCME+MEICA+MGTR)."
  exit 2
fi
Vol2SurfSpec="$VOL2SURF_INPUTS"
if should_run_stage "vol2surf"; then
  run_module "vol2surf" bash "$FUNC_VOL2SURF_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$Vol2SurfSpec" "$START_SESSION" "$AtlasSpace" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
else
  echo "Skipping vol2surf (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Generating Motion & GrayPlots QC"
if should_run_stage "vol2surf"; then
  run_motion_grayplots_qc
  if [[ "$STOP_AFTER_MODULE" == "vol2surf" ]]; then
    echo "Stopping after vol2surf (STOP_AFTER_MODULE=vol2surf)"
    exit 0
  fi
else
  echo "Skipping Motion & GrayPlots QC (START_FROM_MODULE=${START_FROM_MODULE})"
fi

if [[ "${CONCAT_ENABLE}" == "1" ]]; then
  if should_run_stage "concat"; then
    print_section "Running concat module"
    [ -f "$FUNC_CONCAT_MODULE" ] || { echo "ERROR: missing module: $FUNC_CONCAT_MODULE"; exit 2; }
    run_module "concat" bash "$FUNC_CONCAT_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    if [[ "${CRAWLING_SEED_FC_ENABLE}" == "1" ]]; then
      print_section "Generating crawling seed FC movie"
      [ -f "$FUNC_CRAWLING_SEED_FC_MODULE" ] || { echo "ERROR: missing module: $FUNC_CRAWLING_SEED_FC_MODULE"; exit 2; }
      MODULE_DISPLAY_TAG="concat:crawling_seed_fc" MODULE_LOG_TAG="concat_crawling_seed_fc" \
        run_module "crawling_seed_fc" bash "$FUNC_CRAWLING_SEED_FC_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    else
      echo "Skipping crawling seed FC movie (CRAWLING_SEED_FC_ENABLE=${CRAWLING_SEED_FC_ENABLE})"
    fi
    if [[ "$STOP_AFTER_MODULE" == "concat" ]]; then
      echo "Stopping after concat (STOP_AFTER_MODULE=concat)"
      exit 0
    fi
  else
    echo "Skipping concat (START_FROM_MODULE=${START_FROM_MODULE})"
  fi
fi

if [[ "${NSI_ENABLE}" == "1" ]]; then
  if should_run_stage "nsi"; then
    print_section "Running NSI module"
    [ -f "$FUNC_NSI_MODULE" ] || { echo "ERROR: missing module: $FUNC_NSI_MODULE"; exit 2; }
    run_module "nsi" bash "$FUNC_NSI_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    if [[ "$STOP_AFTER_MODULE" == "nsi" ]]; then
      echo "Stopping after nsi (STOP_AFTER_MODULE=nsi)"
      exit 0
    fi
  else
    echo "Skipping nsi (START_FROM_MODULE=${START_FROM_MODULE})"
  fi
fi

if [[ "${PFM_ENABLE}" == "1" ]]; then
  if should_run_stage "pfm"; then
    print_section "Running PFM module"
    [ -f "$FUNC_PFM_MODULE" ] || { echo "ERROR: missing module: $FUNC_PFM_MODULE"; exit 2; }
    run_module "pfm" bash "$FUNC_PFM_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    if [[ "$STOP_AFTER_MODULE" == "pfm" ]]; then
      echo "Stopping after pfm (STOP_AFTER_MODULE=pfm)"
      exit 0
    fi
  else
    echo "Skipping pfm (START_FROM_MODULE=${START_FROM_MODULE})"
  fi
fi

if [[ "${PFM_INFOMAP_UPDATE_ENABLE}" == "1" ]]; then
  if should_run_stage "pfm_update"; then
    print_section "Applying Infomap manual updates"
    [ -f "$FUNC_PFM_INFOMAP_UPDATE_MODULE" ] || { echo "ERROR: missing module: $FUNC_PFM_INFOMAP_UPDATE_MODULE"; exit 2; }
    run_module "pfm_update" bash "$FUNC_PFM_INFOMAP_UPDATE_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    if [[ "$STOP_AFTER_MODULE" == "pfm_update" || "$STOP_AFTER_MODULE" == "infomap_update" ]]; then
      echo "Stopping after pfm_update (STOP_AFTER_MODULE=${STOP_AFTER_MODULE})"
      exit 0
    fi
  else
    echo "Skipping Infomap manual update (START_FROM_MODULE=${START_FROM_MODULE})"
  fi
fi

if [[ "${TANS_ENABLE}" == "1" ]]; then
  if should_run_stage "tans"; then
    print_section "Running TANS module"
    [ -f "$FUNC_TANS_MODULE" ] || { echo "ERROR: missing module: $FUNC_TANS_MODULE"; exit 2; }
    run_module "tans" bash "$FUNC_TANS_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
    if [[ "$STOP_AFTER_MODULE" == "tans" ]]; then
      echo "Stopping after tans (STOP_AFTER_MODULE=tans)"
      exit 0
    fi
  else
    echo "Skipping TANS (START_FROM_MODULE=${START_FROM_MODULE})"
  fi
fi

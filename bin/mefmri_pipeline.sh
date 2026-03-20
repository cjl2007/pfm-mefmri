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
DEFAULT_CONFIG="$SCRIPT_DIR/../config/mefmri_wrapper_config.sh"
ConfigFile="${ConfigFileArg:-${CONFIG_FILE:-$DEFAULT_CONFIG}}"

if [ ! -f "$ConfigFile" ]; then
  echo "ERROR: missing config file: $ConfigFile"
  exit 2
fi

set -a
source "$ConfigFile"
set +a

# Optional config-specified FreeSurfer license path.
if [[ -n "${FS_LICENSE_FILE:-}" ]]; then
  export FS_LICENSE="$FS_LICENSE_FILE"
fi

# Core paths and wrapper controls.
: "${MEDIR:=$(cd "$SCRIPT_DIR/.." && pwd)}"
: "${EnvironmentScript:=$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh}"
: "${START_SESSION:=1}"
: "${START_FROM_MODULE:=validate}" # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|meica|mgtr|vol2surf|concat|nsi|pfm
: "${STOP_AFTER_MODULE:=}" # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|meica|mgtr|vol2surf|concat|nsi|pfm

# Module entrypoints.
: "${VALIDATE_MODULE:=$MEDIR/modules/mefmri_validate_inputs.sh}"
: "${ANAT_HCP_MODULE:=$MEDIR/modules/mefmri_anat_hcp.sh}"
: "${ANAT_CHARM_MODULE:=$MEDIR/modules/mefmri_anat_charm.sh}"
: "${FUNC_FIELDMAPS_MODULE:=$MEDIR/modules/mefmri_func_fieldmaps.sh}"
: "${FUNC_COREG_MODULE:=$MEDIR/modules/mefmri_func_coreg.sh}"
: "${FUNC_HEADMOTION_MODULE:=$MEDIR/modules/mefmri_func_headmotion.sh}"
: "${FUNC_MEICA_MODULE:=$MEDIR/modules/mefmri_func_meica.sh}"
: "${FUNC_MGTR_MODULE:=$MEDIR/modules/mefmri_func_mgtr.sh}"
: "${FUNC_VOL2SURF_MODULE:=$MEDIR/modules/mefmri_func_vol2surf.sh}"
: "${FUNC_CONCAT_MODULE:=$MEDIR/modules/mefmri_func_concat.sh}"
: "${FUNC_NSI_MODULE:=$MEDIR/modules/mefmri_func_nsi.sh}"
: "${FUNC_PFM_MODULE:=$MEDIR/modules/mefmri_func_pfm.sh}"

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
: "${RUN_CONFIG_SNAPSHOT:=1}" # 0|1

# Functional naming/outputs.
: "${VOL2SURF_INPUTS:=}"
: "${FUNC_DIRNAME:=rest}"
: "${FUNC_FILE_PREFIX:=Rest}"

# Module-specific thread controls.
: "${THREADS_DEFAULT:=8}"
: "${THREADS_ANAT_HCP:=$THREADS_DEFAULT}"
: "${THREADS_FIELDMAPS:=$THREADS_DEFAULT}"
: "${THREADS_COREG:=$THREADS_DEFAULT}"
: "${THREADS_HEADMOTION:=$THREADS_DEFAULT}"
: "${THREADS_MEICA:=$THREADS_DEFAULT}"
: "${PIPELINE_QUIET_MODULE_OUTPUT:=1}" # 0|1
: "${PIPELINE_LOG_TAIL_LINES:=40}"

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

stage_index() {
  case "$1" in
    validate) echo 5 ;;
    anat_hcp) echo 10 ;;
    anat_charm) echo 20 ;;
    fieldmaps) echo 30 ;;
    coreg) echo 40 ;;
    headmotion) echo 50 ;;
    meica) echo 60 ;;
    mgtr) echo 70 ;;
    vol2surf) echo 80 ;;
    concat) echo 90 ;;
    nsi) echo 100 ;;
    pfm) echo 110 ;;
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
if should_run_stage "anat_hcp"; then
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
echo "Functional naming: func/${FUNC_DIRNAME}, prefix ${FUNC_FILE_PREFIX}_*"
echo "Threads (anat_hcp,fieldmaps,coreg,headmotion,meica): ${THREADS_ANAT_HCP},${THREADS_FIELDMAPS},${THREADS_COREG},${THREADS_HEADMOTION},${THREADS_MEICA}"
echo "MEICA defaults: tedana_env=${TEDANA_ENV:-unset}, compat=${TEDANA_COMPAT_MODE:-unset}, pca=${MEPCA}"
echo "Masking defaults: CHARM_BRAIN_MASK_MODE=${CHARM_BRAIN_MASK_MODE:-unset}, VOL2SURF_USE_CORTICAL_RIBBON_MASK=${VOL2SURF_USE_CORTICAL_RIBBON_MASK:-unset}"
echo "Reclass defaults: mode=${MEICA_CLASSIFIER_MODE:-unset}, nsi_kill=${MEICA_NSI_KILL_MODE:-unset}, reports_disabled=${MEICA_RECLASS_NO_REPORTS:-unset}"
echo "Module output mode: quiet=${PIPELINE_QUIET_MODULE_OUTPUT}"

MODULE_LOG_DIR="$SubjectDir/func/qa/ModuleLogs"
mkdir -p "$MODULE_LOG_DIR"

run_module() {
  local module="$1"
  shift

  if [[ "${PIPELINE_QUIET_MODULE_OUTPUT}" == "0" ]]; then
    "$@"
    return $?
  fi

  local ts logfile rc
  ts="$(date +%Y%m%d_%H%M%S)"
  logfile="$MODULE_LOG_DIR/${ts}_${module}.log"
  echo "[$module] log: $logfile"

  # Always capture complete output to log; only stream high-level markers to terminal.
  set +e
  "$@" 2>&1 | awk -v logfile_path="$logfile" '
    {
      print >> logfile_path
      fflush(logfile_path)
      if ($0 ~ /^\[(concat|nsi|pfm)\] complete$/) {
        next
      }
      if ($0 ~ /^\[coreg\]/ || $0 ~ /^\[headmotion\]/ || $0 ~ /^\[concat\]/ || $0 ~ /^\[nsi\]/ || $0 ~ /^\[pfm\]/ ||
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
    echo "[$module] FAILED (exit $rc). Log: $logfile" >&2
    echo "[$module] error summary:" >&2
    if ! tail -n 400 "$logfile" | awk '
      /Traceback \(most recent call last\):/ ||
      /^ERROR:/ ||
      /FileNotFoundError:/ ||
      /Exception:/ ||
      /FAILED/ ||
      /can'\''t open file/ { print; found=1 }
      END { exit(found?0:1) }
    ' >&2; then
      echo "[$module] (no explicit traceback/error markers found; showing last ${PIPELINE_LOG_TAIL_LINES} lines)" >&2
      tail -n "$PIPELINE_LOG_TAIL_LINES" "$logfile" >&2 || true
    fi
    return "$rc"
  fi

  echo "[$module] complete"
}

print_section() {
  local title="$1"
  echo
  echo "$title"
}

if [[ "${RUN_CONFIG_SNAPSHOT}" == "1" ]]; then
  RUN_META_DIR="$SubjectDir/func/qa/RunMetadata"
  mkdir -p "$RUN_META_DIR"
  RUN_META_FILE="$RUN_META_DIR/pipeline_run_$(date +%Y%m%d_%H%M%S).txt"
  {
    echo "date=$(date --iso-8601=seconds)"
    echo "subject_dir=$SubjectDir"
    echo "config_file=$ConfigFile"
    echo "start_from_module=$START_FROM_MODULE"
    echo "stop_after_module=${STOP_AFTER_MODULE:-}"
    echo "atlas_space=$AtlasSpace"
    echo "func_dirname=$FUNC_DIRNAME"
    echo "func_file_prefix=$FUNC_FILE_PREFIX"
    echo "charm_brain_mask_mode=${CHARM_BRAIN_MASK_MODE:-}"
    echo "vol2surf_use_cortical_ribbon_mask=${VOL2SURF_USE_CORTICAL_RIBBON_MASK:-}"
    echo "tedana_env=${TEDANA_ENV:-}"
    echo "tedana_compat_mode=${TEDANA_COMPAT_MODE:-}"
    echo "mepca=$MEPCA"
    echo "meica_classifier_mode=${MEICA_CLASSIFIER_MODE:-}"
    echo "meica_nsi_kill_mode=${MEICA_NSI_KILL_MODE:-}"
    echo "meica_reclass_no_reports=${MEICA_RECLASS_NO_REPORTS:-}"
    echo "concat_enable=${CONCAT_ENABLE}"
    echo "nsi_enable=${NSI_ENABLE}"
    echo "pfm_enable=${PFM_ENABLE}"
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
  [ -f "$ANAT_HCP_MODULE" ] || { echo "ERROR: missing module: $ANAT_HCP_MODULE"; exit 2; }
  run_module "anat_hcp" bash "$ANAT_HCP_MODULE" "$StudyFolder" "$Subject" "$THREADS_ANAT_HCP"
  if [[ "$STOP_AFTER_MODULE" == "anat_hcp" ]]; then
    echo "Stopping after anat_hcp (STOP_AFTER_MODULE=anat_hcp)"
    exit 0
  fi
else
  echo "Skipping anat_hcp (START_FROM_MODULE=${START_FROM_MODULE})"
fi
if should_run_stage "anat_charm"; then
  [ -f "$ANAT_CHARM_MODULE" ] || { echo "ERROR: missing module: $ANAT_CHARM_MODULE"; exit 2; }
  if [ -n "$CHARM_BIN" ]; then
    run_module "anat_charm" bash "$ANAT_CHARM_MODULE" "$StudyFolder" "$Subject" "$CHARM_BIN"
  else
    run_module "anat_charm" bash "$ANAT_CHARM_MODULE" "$StudyFolder" "$Subject"
  fi
  if [[ "$STOP_AFTER_MODULE" == "anat_charm" ]]; then
    echo "Stopping after anat_charm (STOP_AFTER_MODULE=anat_charm)"
    exit 0
  fi
else
  echo "Skipping anat_charm (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Processing fieldmaps"
[ -f "$FUNC_FIELDMAPS_MODULE" ] || { echo "ERROR: missing module: $FUNC_FIELDMAPS_MODULE"; exit 2; }
[ -x "$FUNC_FIELDMAPS_MODULE" ] || chmod +x "$FUNC_FIELDMAPS_MODULE"
if should_run_stage "fieldmaps"; then
  run_module "fieldmaps" bash "$FUNC_FIELDMAPS_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$THREADS_FIELDMAPS" "$START_SESSION"
  if [[ "$STOP_AFTER_MODULE" == "fieldmaps" ]]; then
    echo "Stopping after fieldmaps (STOP_AFTER_MODULE=fieldmaps)"
    exit 0
  fi
else
  echo "Skipping fieldmaps (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Coregistering SBrefs to anatomical image"
[ -f "$FUNC_COREG_MODULE" ] || { echo "ERROR: missing module: $FUNC_COREG_MODULE"; exit 2; }
if should_run_stage "coreg"; then
  run_module "coreg" bash "$FUNC_COREG_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$AtlasTemplate" "$DOF" "$THREADS_COREG" "$START_SESSION" "$AtlasSpace" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
  if [[ "$STOP_AFTER_MODULE" == "coreg" ]]; then
    echo "Stopping after coreg (STOP_AFTER_MODULE=coreg)"
    exit 0
  fi
else
  echo "Skipping coreg (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Applying slice-time/headmotion/distortion corrections"
[ -f "$FUNC_HEADMOTION_MODULE" ] || { echo "ERROR: missing module: $FUNC_HEADMOTION_MODULE"; exit 2; }
if should_run_stage "headmotion"; then
  run_module "headmotion" bash "$FUNC_HEADMOTION_MODULE" "$MEDIR" "$Subject" "$StudyFolder" "$AtlasTemplate" "$DOF" "$THREADS_HEADMOTION" "$START_SESSION" "$AtlasSpace" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
  if [[ "$STOP_AFTER_MODULE" == "headmotion" ]]; then
    echo "Stopping after headmotion (STOP_AFTER_MODULE=headmotion)"
    exit 0
  fi
else
  echo "Skipping headmotion (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Running ME-ICA denoising"
[ -f "$FUNC_MEICA_MODULE" ] || { echo "ERROR: missing module: $FUNC_MEICA_MODULE"; exit 2; }
if should_run_stage "meica"; then
  run_module "meica" bash "$FUNC_MEICA_MODULE" "$Subject" "$StudyFolder" "$THREADS_MEICA" "$MEPCA" "$MaxIterations" "$MaxRestarts" "$START_SESSION" "$MEDIR"
  if [[ "$STOP_AFTER_MODULE" == "meica" ]]; then
    echo "Stopping after meica (STOP_AFTER_MODULE=meica)"
    exit 0
  fi
else
  echo "Skipping meica (START_FROM_MODULE=${START_FROM_MODULE})"
fi

print_section "Running MGTR"
[ -f "$FUNC_MGTR_MODULE" ] || { echo "ERROR: missing module: $FUNC_MGTR_MODULE"; exit 2; }
if should_run_stage "mgtr"; then
  run_module "mgtr" bash "$FUNC_MGTR_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION"
  if [[ "$STOP_AFTER_MODULE" == "mgtr" ]]; then
    echo "Stopping after mgtr (STOP_AFTER_MODULE=mgtr)"
    exit 0
  fi
else
  echo "Skipping mgtr (START_FROM_MODULE=${START_FROM_MODULE})"
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
  if [[ "$STOP_AFTER_MODULE" == "vol2surf" ]]; then
    echo "Stopping after vol2surf (STOP_AFTER_MODULE=vol2surf)"
    exit 0
  fi
else
  echo "Skipping vol2surf (START_FROM_MODULE=${START_FROM_MODULE})"
fi

if [[ "${CONCAT_ENABLE}" == "1" ]]; then
  if should_run_stage "concat"; then
    print_section "Running concat module"
    [ -f "$FUNC_CONCAT_MODULE" ] || { echo "ERROR: missing module: $FUNC_CONCAT_MODULE"; exit 2; }
    run_module "concat" bash "$FUNC_CONCAT_MODULE" "$Subject" "$StudyFolder" "$MEDIR" "$START_SESSION" "$FUNC_DIRNAME" "$FUNC_FILE_PREFIX"
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

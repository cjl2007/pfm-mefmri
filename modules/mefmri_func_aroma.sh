#!/usr/bin/env bash
# Classic ICA-AROMA denoising.
#
# A future experimental NSI/surface reclassification path is retained behind
# AROMA_NSI_RECLASSIFY_ENABLE=1, but the default behavior is classic AROMA.

set -euo pipefail
IFS=$'\n\t'

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
NTHREADS="${3:?missing NTHREADS}"
StartSession="${4:?missing StartSession}"
MEDIR="${5:?missing MEDIR}"

Subdir="$StudyFolder/$Subject"
FuncDirName="${FUNC_DIRNAME:-rest}"
FuncFilePrefix="${FUNC_FILE_PREFIX:-Rest}"
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-$FuncDirName}"
AtlasSpace="${AtlasSpace:-T1w}"

AROMA_ENV="${AROMA_ENV:-aroma}"
AROMA_ACTIVATE_MODE="${AROMA_ACTIVATE_MODE:-conda_activate}" # conda_activate|conda_run|direct
AROMA_PARALLEL_JOBS="${AROMA_PARALLEL_JOBS:-$NTHREADS}"
AROMA_OUT_SUBDIR="${AROMA_OUT_SUBDIR:-Aroma}"
AROMA_OVERWRITE="${AROMA_OVERWRITE:-1}"
AROMA_PLOT_REPORTS="${AROMA_PLOT_REPORTS:-0}"
AROMA_FEATURES_OUTPUT="${AROMA_FEATURES_OUTPUT:-1}"
AROMA_CLEAN_TYPE="${AROMA_CLEAN_TYPE:-nonaggr}" # nonaggr|aggr
AROMA_BIN_OVERRIDE="${AROMA_BIN_OVERRIDE:-}"
AROMA_PYTHON="${AROMA_PYTHON:-python}"
AROMA_SCREEN_PYTHON="${AROMA_SCREEN_PYTHON:-${PIPELINE_PYTHON:-python3}}"
AROMA_NSI_RECLASSIFY_ENABLE="${AROMA_NSI_RECLASSIFY_ENABLE:-0}"
AROMA_NSI_THRESHOLD="${AROMA_NSI_THRESHOLD:-0.05}"
SINGLE_ECHO_ECHO_INDEX="${SINGLE_ECHO_ECHO_INDEX:-1}"
PIPELINE_SOURCE_FUNC_TAG="${PIPELINE_SOURCE_FUNC_TAG:-E${SINGLE_ECHO_ECHO_INDEX}}"
PIPELINE_DENOISE_OUTPUT_TAG="${PIPELINE_DENOISE_OUTPUT_TAG:-${PIPELINE_SOURCE_FUNC_TAG}+AROMA}"
NETWORK_PRIORS_MAT="${NETWORK_PRIORS_MAT:-$MEDIR/res0urces/Priors.mat}"
AROMA_SCREEN_PY="${MEDIR}/lib/aroma_screen_components.py"

log() { echo "[aroma] $*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

case "$AROMA_ACTIVATE_MODE" in
  conda_activate|conda_run|direct) ;;
  *) die "Invalid AROMA_ACTIVATE_MODE=$AROMA_ACTIVATE_MODE (use conda_activate|conda_run|direct)" ;;
esac
case "$AROMA_OVERWRITE" in
  0|1) ;;
  *) die "Invalid AROMA_OVERWRITE=$AROMA_OVERWRITE (use 0|1)" ;;
esac
case "$AROMA_PLOT_REPORTS" in
  0|1) ;;
  *) die "Invalid AROMA_PLOT_REPORTS=$AROMA_PLOT_REPORTS (use 0|1)" ;;
esac
case "$AROMA_FEATURES_OUTPUT" in
  0|1) ;;
  *) die "Invalid AROMA_FEATURES_OUTPUT=$AROMA_FEATURES_OUTPUT (use 0|1)" ;;
esac
case "$AROMA_CLEAN_TYPE" in
  nonaggr|aggr) ;;
  *) die "Invalid AROMA_CLEAN_TYPE=$AROMA_CLEAN_TYPE (use nonaggr|aggr)" ;;
esac
case "$AROMA_NSI_RECLASSIFY_ENABLE" in
  0|1) ;;
  *) die "Invalid AROMA_NSI_RECLASSIFY_ENABLE=$AROMA_NSI_RECLASSIFY_ENABLE (use 0|1)" ;;
esac

need_cmd fsl_regfilt
if [[ "$AROMA_NSI_RECLASSIFY_ENABLE" == "1" ]]; then
  need_cmd wb_command
fi
if [[ "$AROMA_ACTIVATE_MODE" != "direct" ]]; then
  need_cmd conda
fi

resolve_aroma_bin() {
  if [[ -n "$AROMA_BIN_OVERRIDE" ]]; then
    [[ -f "$AROMA_BIN_OVERRIDE" ]] || die "AROMA_BIN_OVERRIDE does not exist: $AROMA_BIN_OVERRIDE"
    echo "$AROMA_BIN_OVERRIDE"
    return 0
  fi
  if [[ "$AROMA_ACTIVATE_MODE" == "direct" ]]; then
    command -v ICA_AROMA.py || true
    return 0
  fi
  local base
  base="$(conda info --base 2>/dev/null || true)"
  [[ -n "${base:-}" ]] || { echo ""; return 0; }
  local env_dir="$base/envs/$AROMA_ENV"
  [[ -d "$env_dir" ]] || { echo ""; return 0; }

  # Prefer resolving via the activated environment PATH (covers pip installs that
  # add ICA_AROMA.py as an entrypoint script).
  local aroma_path=""
  aroma_path="$(conda run -n "$AROMA_ENV" sh -c 'command -v ICA_AROMA.py || true' 2>/dev/null || true)"
  if [[ -n "$aroma_path" && -f "$aroma_path" ]]; then
    echo "$aroma_path"
    return 0
  fi

  # Fallback: common on-disk location.
  if [[ -x "$env_dir/bin/ICA_AROMA.py" ]]; then
    echo "$env_dir/bin/ICA_AROMA.py"
    return 0
  fi
  if [[ -f "$env_dir/bin/ICA_AROMA.py" ]]; then
    echo "$env_dir/bin/ICA_AROMA.py"
    return 0
  fi
  echo ""
}

AROMA_BIN="$(resolve_aroma_bin)"
if [[ -z "$AROMA_BIN" ]]; then
  base="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${base:-}" && ! -d "$base/envs/$AROMA_ENV" ]]; then
    die "ICA-AROMA conda env not found: $AROMA_ENV (looked in $base/envs). Set AROMA_ENV to an existing env or set AROMA_BIN_OVERRIDE to ICA_AROMA.py."
  fi
  die "Could not resolve ICA_AROMA.py for env=$AROMA_ENV mode=$AROMA_ACTIVATE_MODE. Ensure ICA-AROMA is installed so that 'ICA_AROMA.py' is on PATH in that env, or set AROMA_BIN_OVERRIDE."
fi
if [[ "$AROMA_NSI_RECLASSIFY_ENABLE" == "1" ]]; then
  [[ -f "$AROMA_SCREEN_PY" ]] || die "Missing aroma screening backend: $AROMA_SCREEN_PY"
  [[ -f "$NETWORK_PRIORS_MAT" ]] || die "Missing Priors.mat for ICA-AROMA NSI screening: $NETWORK_PRIORS_MAT"
fi

aroma_exec() {
  local bin="$1"
  shift
  case "$AROMA_ACTIVATE_MODE" in
    conda_run)
      conda run -n "$AROMA_ENV" "$AROMA_PYTHON" "$bin" "$@"
      ;;
    conda_activate)
      bash -lc '
        set -euo pipefail
        BASE="$(conda info --base 2>/dev/null || true)"
        [[ -n "${BASE:-}" && -f "$BASE/etc/profile.d/conda.sh" ]] && source "$BASE/etc/profile.d/conda.sh"
        conda activate "'"$AROMA_ENV"'"
        exec "'"$AROMA_PYTHON"'" "'"$bin"'" "$@"
      ' _ "$@"
      ;;
    direct)
      "$AROMA_PYTHON" "$bin" "$@"
      ;;
  esac
}

map_qc_volume_to_cifti() {
  local in_vol="$1"
  local out_cifti="$2"
  local tmpdir="$3"
  local dim4="$4"

  local pial_l="$Subdir/anat/T1w/Native/${Subject}.L.pial.native.surf.gii"
  local white_l="$Subdir/anat/T1w/Native/${Subject}.L.white.native.surf.gii"
  local mid_l="$Subdir/anat/T1w/Native/${Subject}.L.midthickness.native.surf.gii"
  local pial_r="$Subdir/anat/T1w/Native/${Subject}.R.pial.native.surf.gii"
  local white_r="$Subdir/anat/T1w/Native/${Subject}.R.white.native.surf.gii"
  local mid_r="$Subdir/anat/T1w/Native/${Subject}.R.midthickness.native.surf.gii"
  local roi_l="$Subdir/anat/MNINonLinear/Native/${Subject}.L.roi.native.shape.gii"
  local roi_r="$Subdir/anat/MNINonLinear/Native/${Subject}.R.roi.native.shape.gii"
  local roi32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.atlasroi.32k_fs_LR.shape.gii"
  local roi32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.atlasroi.32k_fs_LR.shape.gii"
  local reg_l="$Subdir/anat/MNINonLinear/Native/${Subject}.L.sphere.MSMSulc.native.surf.gii"
  local reg_r="$Subdir/anat/MNINonLinear/Native/${Subject}.R.sphere.MSMSulc.native.surf.gii"
  local reg32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.sphere.32k_fs_LR.surf.gii"
  local reg32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.sphere.32k_fs_LR.surf.gii"
  local mid32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.midthickness.32k_fs_LR.surf.gii"
  local mid32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.midthickness.32k_fs_LR.surf.gii"
  local subcort="$Subdir/func/rois/Subcortical_ROIs_acpc.nii.gz"

  rm -rf "$tmpdir" 2>/dev/null || true
  mkdir -p "$tmpdir"
  local l_native="$tmpdir/lh.native.shape.gii"
  local r_native="$tmpdir/rh.native.shape.gii"
  local l_32k="$tmpdir/lh.32k_fs_LR.shape.gii"
  local r_32k="$tmpdir/rh.32k_fs_LR.shape.gii"

  wb_command -volume-to-surface-mapping "$in_vol" "$mid_l" "$l_native" -ribbon-constrained "$white_l" "$pial_l"
  wb_command -volume-to-surface-mapping "$in_vol" "$mid_r" "$r_native" -ribbon-constrained "$white_r" "$pial_r"
  wb_command -metric-dilate "$l_native" "$mid_l" 10 "$l_native" -nearest
  wb_command -metric-dilate "$r_native" "$mid_r" 10 "$r_native" -nearest
  wb_command -metric-mask "$l_native" "$roi_l" "$l_native"
  wb_command -metric-mask "$r_native" "$roi_r" "$r_native"
  wb_command -metric-resample "$l_native" "$reg_l" "$reg32_l" ADAP_BARY_AREA "$l_32k" -area-surfs "$mid_l" "$mid32_l" -current-roi "$roi_l"
  wb_command -metric-resample "$r_native" "$reg_r" "$reg32_r" ADAP_BARY_AREA "$r_32k" -area-surfs "$mid_r" "$mid32_r" -current-roi "$roi_r"
  wb_command -metric-mask "$l_32k" "$roi32_l" "$l_32k"
  wb_command -metric-mask "$r_32k" "$roi32_r" "$r_32k"

  if [[ "$dim4" -gt 1 ]]; then
    wb_command -cifti-create-dense-timeseries "$out_cifti" \
      -volume "$in_vol" "$subcort" \
      -left-metric "$l_32k" -roi-left "$roi32_l" \
      -right-metric "$r_32k" -roi-right "$roi32_r"
  else
    wb_command -cifti-create-dense-scalar "$out_cifti" \
      -volume "$in_vol" "$subcort" \
      -left-metric "$l_32k" \
      -right-metric "$r_32k"
  fi
  rm -rf "$tmpdir" 2>/dev/null || true
}

component_filter_csv() {
  local path="$1"
  [[ -f "$path" ]] || return 0
  tr ', \n\t' ',' < "$path" | sed 's/,,*/,/g; s/^,//; s/,$//'
}

list_runs() {
  find -L "$Subdir/func/$FuncDirName" -mindepth 2 -maxdepth 2 -type d -name 'run_*' | sort -V
}

process_run() {
  local run_dir="$1"
  local rel="${run_dir#$Subdir/func/$FuncDirName/}"
  local input_vol="$run_dir/${FuncFilePrefix}_${PIPELINE_SOURCE_FUNC_TAG}.nii.gz"
  local input_vol_acpc="$run_dir/${FuncFilePrefix}_${PIPELINE_SOURCE_FUNC_TAG}_acpc.nii.gz"
  local out_dir="$run_dir/${AROMA_OUT_SUBDIR}"
  local mask=""
  if [[ "$AtlasSpace" == "MNINonlinear" ]]; then
    mask="$Subdir/func/xfms/$FuncXfmsDir/T1w_nonlin_brain_func_mask.nii.gz"
  else
    mask="$Subdir/func/xfms/$FuncXfmsDir/T1w_acpc_brain_func_mask.nii.gz"
  fi
  local warp="$Subdir/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz"
  local motion="$run_dir/MCF.par"
  local tr_file="$run_dir/TR.txt"
  local classified_motion="$out_dir/classified_motion_ICs.txt"
  local comp_nii="$out_dir/melodic.ica/melodic_IC.nii.gz"
  local mix_tsv="$out_dir/melodic.ica/melodic_mix"
  local comp_cifti="$out_dir/melodic_IC.dtseries.nii"
  local screen_dir="$out_dir/ComponentScreen"
  local accepted_1b="$screen_dir/AcceptedComponents_1based.txt"
  local rejected_1b="$screen_dir/RejectedComponents_1based.txt"
  local final_reject_csv=""
  local denoised_out="$run_dir/${FuncFilePrefix}_${PIPELINE_DENOISE_OUTPUT_TAG}.nii.gz"
  local tmpdir="$run_dir/.tmp_aroma_cifti"

  if [[ ! -f "$input_vol" && -f "$input_vol_acpc" ]]; then
    cp -f "$input_vol_acpc" "$input_vol"
  fi
  [[ -f "$input_vol" ]] || die "ICA-AROMA input missing for $rel: $input_vol"
  [[ -f "$tr_file" ]] || die "ICA-AROMA TR missing for $rel: $tr_file"
  [[ -f "$motion" ]] || {
    motion="$run_dir/mcf.par"
    [[ -f "$motion" ]] || die "ICA-AROMA motion file missing for $rel: $run_dir/MCF.par"
  }
  if [[ ! -f "$mask" && -f "$Subdir/func/xfms/$FuncXfmsDir/T1w_func_brain_mask.nii.gz" ]]; then
    mask="$Subdir/func/xfms/$FuncXfmsDir/T1w_func_brain_mask.nii.gz"
  fi
  [[ -f "$mask" ]] || die "ICA-AROMA brain mask missing: $mask"

  rm -rf "$out_dir"
  mkdir -p "$out_dir"
  local tr
  tr="$(<"$tr_file")"

  log "start $rel"
  local -a aroma_args=(
    -i "$input_vol"
    -o "$out_dir"
    -mc "$motion"
    -tr "$tr"
    -m "$mask"
  )
  if [[ "$AtlasSpace" == "T1w" && -f "$warp" ]]; then
    aroma_args+=( -w "$warp" )
  fi
  if [[ "$AROMA_OVERWRITE" == "1" ]]; then
    aroma_args+=( -overwrite )
  fi

  aroma_exec "$AROMA_BIN" "${aroma_args[@]}"

  [[ -f "$classified_motion" ]] || die "ICA-AROMA did not produce classified_motion_ICs.txt for $rel"
  [[ -f "$comp_nii" ]] || die "ICA-AROMA did not produce melodic_IC.nii.gz for $rel"
  [[ -f "$mix_tsv" ]] || die "ICA-AROMA did not produce melodic_mix for $rel"

  if [[ "$AROMA_NSI_RECLASSIFY_ENABLE" == "1" ]]; then
    map_qc_volume_to_cifti "$comp_nii" "$comp_cifti" "$tmpdir" "$(fslval "$comp_nii" dim4)"

    "$AROMA_SCREEN_PYTHON" "$AROMA_SCREEN_PY" \
      --aroma-dir "$out_dir" \
      --betas-cifti "$comp_cifti" \
      --priors-mat "$NETWORK_PRIORS_MAT" \
      --nsi-threshold "$AROMA_NSI_THRESHOLD" \
      --kill-priority-enable "${AROMA_NSI_KILL_PRIORITY_ENABLE:-1}" \
      --kill-priority-w-nsi "${AROMA_NSI_KILL_PRIORITY_W_NSI:-0.60}" \
      --kill-priority-w-var "${AROMA_NSI_KILL_PRIORITY_W_VAR:-0.25}" \
      --kill-priority-w-aroma "${AROMA_NSI_KILL_PRIORITY_W_AROMA:-0.15}" \
      --kill-var-floor-quantile "${AROMA_NSI_KILL_VAR_FLOOR_QUANTILE:-0.60}" \
      --kill-cumvar-cap "${AROMA_NSI_KILL_CUMVAR_CAP:-0.95}" \
      --kill-max-frac "${AROMA_NSI_KILL_MAX_FRAC:-1.00}" \
      --kill-max-count "${AROMA_NSI_KILL_MAX_COUNT:-0}" \
      --motion-priority-enable "${AROMA_MOTION_PRIORITY_ENABLE:-0}" \
      --motion-remove-frac "${AROMA_MOTION_REMOVE_FRAC:-1.00}" \
      --motion-var-floor-quantile "${AROMA_MOTION_VAR_FLOOR_QUANTILE:-0.60}" \
      --motion-cumvar-cap "${AROMA_MOTION_CUMVAR_CAP:-1.00}" \
      --motion-priority-w-var "${AROMA_MOTION_PRIORITY_W_VAR:-0.50}" \
      --motion-priority-w-aroma "${AROMA_MOTION_PRIORITY_W_AROMA:-0.50}" \
      --make-plot "${AROMA_SCREEN_MAKE_PLOT:-1}" \
      --out-dir "$screen_dir"
    if [[ -s "$rejected_1b" ]]; then
      final_reject_csv="$(paste -sd, "$rejected_1b")"
    fi
  else
    mkdir -p "$screen_dir"
    {
      echo "Classic ICA-AROMA mode."
      echo "AROMA_NSI_RECLASSIFY_ENABLE=0, so NSI/surface component reclassification was not applied."
      echo "Future experimental AROMA reclassification can be enabled explicitly after validation."
    } > "$screen_dir/README.txt"
    final_reject_csv="$(component_filter_csv "$classified_motion")"
  fi

  if [[ -n "$final_reject_csv" ]]; then
    if [[ "$AROMA_CLEAN_TYPE" == "aggr" ]]; then
      fsl_regfilt --in="$input_vol" --design="$mix_tsv" --filter="$final_reject_csv" --out="$denoised_out" -a
    else
      fsl_regfilt --in="$input_vol" --design="$mix_tsv" --filter="$final_reject_csv" --out="$denoised_out"
    fi
  else
    cp -f "$input_vol" "$denoised_out"
  fi

  log "done $rel"
}

log "start subject=${Subject} env=${AROMA_ENV} mode=${AROMA_ACTIVATE_MODE} method=${AROMA_CLEAN_TYPE} nsi_reclassify=${AROMA_NSI_RECLASSIFY_ENABLE}"

mapfile -t RUNS < <(list_runs)
[[ "${#RUNS[@]}" -gt 0 ]] || die "No run directories found in $Subdir/func/$FuncDirName"

for r in "${RUNS[@]}"; do
  session_num="${r%/*}"
  session_num="${session_num##*/}"
  session_num="${session_num#session_}"
  (( session_num >= StartSession )) || continue
  process_run "$r"
done

log "all done subject=${Subject}"

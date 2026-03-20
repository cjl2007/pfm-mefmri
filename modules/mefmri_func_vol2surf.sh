#!/usr/bin/env bash
# CJL; (cjl2007@med.cornell.edu)
# Volume-to-surface mapping for selected functional outputs.

set -euo pipefail
IFS=$'\n\t'

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
Subdir="$StudyFolder/$Subject"
MEDIR="${3:?missing MEDIR}"
InputSpec="${4:?missing InputSpec}"
StartSession="${5:?missing StartSession}"
AtlasSpace="${6:-${AtlasSpace:-T1w}}"
FuncDirName="${7:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${8:-${FUNC_FILE_PREFIX:-Rest}}"
VOL2SURF_USE_CORTICAL_RIBBON_MASK="${VOL2SURF_USE_CORTICAL_RIBBON_MASK:-1}"
VOL2SURF_CIFTI_STAMP="${VOL2SURF_CIFTI_STAMP:-}"

# Coreg currently writes transform-space masks under func/xfms/rest.
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-rest}"

log() { echo "[vol2surf] $*"; }
die() { echo "ERROR: $*" >&2; exit 2; }

case "${AtlasSpace}" in
  T1w|MNINonlinear) ;;
  *) die "invalid AtlasSpace='${AtlasSpace}' (expected T1w or MNINonlinear)" ;;
esac
case "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" in
  0|1) ;;
  *) die "VOL2SURF_USE_CORTICAL_RIBBON_MASK must be 0 or 1 (got '$VOL2SURF_USE_CORTICAL_RIBBON_MASK')" ;;
esac
[[ -n "$InputSpec" ]] || die "missing vol2surf input spec (comma-separated list expected)"

log "AtlasSpace=${AtlasSpace}"
log "Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"
log "Use cortical ribbon mask: ${VOL2SURF_USE_CORTICAL_RIBBON_MASK}"
if [[ -n "$VOL2SURF_CIFTI_STAMP" ]]; then
  log "CIFTI output stamp: ${VOL2SURF_CIFTI_STAMP}"
fi

ensure_cortical_ribbon_func_mask() {
  local subdir="$1"
  local ref_vol="$2"
  local med="$3"
  local atlas_space="$4"
  local out_mask

  if [[ "$atlas_space" == "MNINonlinear" ]]; then
    out_mask="$subdir/func/xfms/$FuncXfmsDir/CorticalRibbon_nonlin_func_mask.nii.gz"
  else
    out_mask="$subdir/func/xfms/$FuncXfmsDir/CorticalRibbon_acpc_func_mask.nii.gz"
  fi

  if [[ -f "$out_mask" ]]; then
    echo "$out_mask"
    return 0
  fi

  local src_mask="$subdir/anat/T1w/CorticalRibbon.nii.gz"
  if [[ ! -f "$src_mask" && -f "$subdir/anat/T1w/CorticalRibbon.ni.gz" ]]; then
    src_mask="$subdir/anat/T1w/CorticalRibbon.ni.gz"
  fi
  if [[ ! -f "$src_mask" ]]; then
    log "WARN: missing CorticalRibbon source; using T1w_acpc_brain_mask fallback"
    src_mask="$subdir/anat/T1w/T1w_acpc_brain_mask.nii.gz"
    [[ -f "$src_mask" ]] || die "missing fallback mask: $src_mask"
  fi

  if [[ "$atlas_space" == "MNINonlinear" ]]; then
    local nonlin_warp="$subdir/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz"
    [[ -f "$nonlin_warp" ]] || die "missing nonlinear warp for AtlasSpace=MNINonlinear: $nonlin_warp"
    applywarp --interp=nn --in="$src_mask" --ref="$ref_vol" --warp="$nonlin_warp" --out="$out_mask"
  else
    flirt -interp nearestneighbour -in "$src_mask" -ref "$ref_vol" -out "$out_mask" -applyxfm -init "$med/res0urces/ident.mat"
  fi

  fslmaths "$out_mask" -bin "$out_mask"
  echo "$out_mask"
}

# Parse comma-separated VOL2SURF_INPUTS.
mapfile -t CIFTI_TAGS < <(echo "$InputSpec" | tr ',' '\n' | sed 's/^ *//; s/ *$//' | awk 'NF')
[[ "${#CIFTI_TAGS[@]}" -gt 0 ]] || die "no valid tags parsed from InputSpec: '$InputSpec'"

mapfile -t SESSION_DIRS < <(
  find "$Subdir/func/$FuncDirName" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V
)

for SES_DIR in "${SESSION_DIRS[@]}"; do
  s="${SES_DIR##*/}"
  s="${s#session_}"
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  (( s >= StartSession )) || continue

  mapfile -t RUN_DIRS < <(
    find "$SES_DIR" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V
  )

  for OUT_DIR in "${RUN_DIRS[@]}"; do
    r="${OUT_DIR##*/}"
    r="${r#run_}"
    [[ "$r" =~ ^[0-9]+$ ]] || continue

    for c in "${CIFTI_TAGS[@]}"; do
      tag="$c"
      tag="${tag#REST_}"
      tag="${tag#Rest_}"

      in_nii="$OUT_DIR/${FuncFilePrefix}_${tag}.nii.gz"
      [[ -f "$in_nii" ]] || { log "WARNING: missing vol2surf input, skipping tag '$c' ($in_nii)"; continue; }

      out_suffix=""
      if [[ -n "$VOL2SURF_CIFTI_STAMP" ]]; then
        out_suffix="_$VOL2SURF_CIFTI_STAMP"
      fi
      out_dt="$OUT_DIR/${FuncFilePrefix}_${tag}${out_suffix}.dtseries.nii"

      CorticalRibbonMask=""
      if [[ "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" == "1" ]]; then
        CorticalRibbonMask="$(ensure_cortical_ribbon_func_mask "$Subdir" "$in_nii" "$MEDIR" "$AtlasSpace")"
      fi

      for hemisphere in lh rh; do
        if [[ "$hemisphere" == "lh" ]]; then
          Hemisphere="L"
        else
          Hemisphere="R"
        fi

        PIAL="$Subdir/anat/T1w/Native/$Subject.$Hemisphere.pial.native.surf.gii"
        WHITE="$Subdir/anat/T1w/Native/$Subject.$Hemisphere.white.native.surf.gii"
        MIDTHICK="$Subdir/anat/T1w/Native/$Subject.$Hemisphere.midthickness.native.surf.gii"
        MIDTHICK_FSLR32k="$Subdir/anat/T1w/fsaverage_LR32k/$Subject.$Hemisphere.midthickness.32k_fs_LR.surf.gii"
        ROI="$Subdir/anat/MNINonLinear/Native/$Subject.$Hemisphere.roi.native.shape.gii"
        ROI_FSLR32k="$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.$Hemisphere.atlasroi.32k_fs_LR.shape.gii"
        REG_MSMSulc="$Subdir/anat/MNINonLinear/Native/$Subject.$Hemisphere.sphere.MSMSulc.native.surf.gii"
        REG_MSMSulc_FSLR32k="$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.$Hemisphere.sphere.32k_fs_LR.surf.gii"

        if [[ "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" == "1" ]]; then
          wb_command -volume-to-surface-mapping "$in_nii" "$MIDTHICK" "$OUT_DIR/$hemisphere.native.shape.gii" \
            -ribbon-constrained "$WHITE" "$PIAL" -volume-roi "$CorticalRibbonMask"
        else
          wb_command -volume-to-surface-mapping "$in_nii" "$MIDTHICK" "$OUT_DIR/$hemisphere.native.shape.gii" \
            -ribbon-constrained "$WHITE" "$PIAL"
        fi

        wb_command -metric-dilate "$OUT_DIR/$hemisphere.native.shape.gii" "$MIDTHICK" 10 "$OUT_DIR/$hemisphere.native.shape.gii" -nearest
        wb_command -metric-mask "$OUT_DIR/$hemisphere.native.shape.gii" "$ROI" "$OUT_DIR/$hemisphere.native.shape.gii"
        wb_command -metric-resample "$OUT_DIR/$hemisphere.native.shape.gii" "$REG_MSMSulc" "$REG_MSMSulc_FSLR32k" \
          ADAP_BARY_AREA "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" -area-surfs "$MIDTHICK" "$MIDTHICK_FSLR32k" -current-roi "$ROI"
        wb_command -metric-mask "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" "$ROI_FSLR32k" "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii"
      done

      tr=$(cat "$OUT_DIR/TR.txt")
      wb_command -cifti-create-dense-timeseries "$out_dt" -volume "$in_nii" "$Subdir/func/rois/Subcortical_ROIs_acpc.nii.gz" \
        -left-metric "$OUT_DIR/lh.32k_fs_LR.shape.gii" -roi-left "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.L.atlasroi.32k_fs_LR.shape.gii" \
        -right-metric "$OUT_DIR/rh.32k_fs_LR.shape.gii" -roi-right "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.R.atlasroi.32k_fs_LR.shape.gii" \
        -timestep "$tr"

      rm -f "$OUT_DIR"/*shape*
    done
  done
done

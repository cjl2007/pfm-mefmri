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
VOL2SURF_USE_GOOD_VOXELS_MASK="${VOL2SURF_USE_GOOD_VOXELS_MASK:-${VOL2SURF_USE_GOOD_VOXELS_MASK_FINAL:-0}}"
VOL2SURF_GOOD_VOXELS_FACTOR="${VOL2SURF_GOOD_VOXELS_FACTOR:-0.5}"
VOL2SURF_GOOD_VOXELS_SIGMA_MM="${VOL2SURF_GOOD_VOXELS_SIGMA_MM:-5}"
VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES="${VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES:-0}"

# Coreg writes transform-space masks under func/xfms/$FuncXfmsDir.
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-$FuncDirName}"

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
case "$VOL2SURF_USE_GOOD_VOXELS_MASK" in
  0|1) ;;
  *) die "VOL2SURF_USE_GOOD_VOXELS_MASK must be 0 or 1 (got '$VOL2SURF_USE_GOOD_VOXELS_MASK')" ;;
esac
case "$VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES" in
  0|1) ;;
  *) die "VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES must be 0 or 1 (got '$VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES')" ;;
esac
[[ -n "$InputSpec" ]] || die "missing vol2surf input spec (comma-separated list expected)"

log "AtlasSpace=${AtlasSpace}"
log "Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"
log "Use cortical ribbon mask: ${VOL2SURF_USE_CORTICAL_RIBBON_MASK}"
log "Use GoodVoxels mask during vol2surf: ${VOL2SURF_USE_GOOD_VOXELS_MASK}"
if [[ "$VOL2SURF_USE_GOOD_VOXELS_MASK" == "1" ]]; then
  log "GoodVoxels defaults: factor=${VOL2SURF_GOOD_VOXELS_FACTOR}, sigma_mm=${VOL2SURF_GOOD_VOXELS_SIGMA_MM}, keep_intermediates=${VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES}"
fi
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

cleanup_goodvoxels_workspace() {
  local gv_dir="$1"
  local goodvoxels_mask="$2"
  local cortical_sampling_mask="$3"
  local subcortical_sampling_mask="$4"

  [[ "$VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES" == "1" ]] && return 0

  find "$gv_dir" -maxdepth 1 -type f ! -name "$(basename "$goodvoxels_mask")" \
    ! -name "$(basename "$cortical_sampling_mask")" \
    ! -name "$(basename "$subcortical_sampling_mask")" \
    ! -name '*cov_before_mask*.shape.gii' \
    ! -name '*cov_after_mask*.shape.gii' \
    ! -name '*cov_before_mask*.dscalar.nii' \
    ! -name '*cov_after_mask*.dscalar.nii' -delete
}

write_goodvoxels_surface_qc() {
  local gv_dir="$1"
  local cov_vol="$2"
  local goodvoxels_mask="$3"
  local subcortical_rois="$4"
  local subcortical_sampling_mask="$5"
  local subdir="$6"
  local subject="$7"

  local hemisphere Hemisphere pial white midthick midthick_fsLR32k roi roi_fsLR32k reg_native reg_fsLR32k
  for hemisphere in lh rh; do
    if [[ "$hemisphere" == "lh" ]]; then
      Hemisphere="L"
    else
      Hemisphere="R"
    fi

    pial="$subdir/anat/T1w/Native/$subject.$Hemisphere.pial.native.surf.gii"
    white="$subdir/anat/T1w/Native/$subject.$Hemisphere.white.native.surf.gii"
    midthick="$subdir/anat/T1w/Native/$subject.$Hemisphere.midthickness.native.surf.gii"
    midthick_fsLR32k="$subdir/anat/T1w/fsaverage_LR32k/$subject.$Hemisphere.midthickness.32k_fs_LR.surf.gii"
    roi="$subdir/anat/MNINonLinear/Native/$subject.$Hemisphere.roi.native.shape.gii"
    roi_fsLR32k="$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.$Hemisphere.atlasroi.32k_fs_LR.shape.gii"
    reg_native="$subdir/anat/MNINonLinear/Native/$subject.$Hemisphere.sphere.MSMSulc.native.surf.gii"
    reg_fsLR32k="$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.$Hemisphere.sphere.32k_fs_LR.surf.gii"

    wb_command -volume-to-surface-mapping "$cov_vol" "$midthick" "$gv_dir/$hemisphere.cov_before_mask.native.shape.gii" \
      -ribbon-constrained "$white" "$pial"
    wb_command -metric-mask "$gv_dir/$hemisphere.cov_before_mask.native.shape.gii" "$roi" "$gv_dir/$hemisphere.cov_before_mask.native.shape.gii"
    wb_command -metric-resample "$gv_dir/$hemisphere.cov_before_mask.native.shape.gii" "$reg_native" "$reg_fsLR32k" \
      ADAP_BARY_AREA "$gv_dir/$hemisphere.cov_before_mask.32k_fs_LR.shape.gii" -area-surfs "$midthick" "$midthick_fsLR32k" -current-roi "$roi"
    wb_command -metric-mask "$gv_dir/$hemisphere.cov_before_mask.32k_fs_LR.shape.gii" "$roi_fsLR32k" "$gv_dir/$hemisphere.cov_before_mask.32k_fs_LR.shape.gii"

    wb_command -volume-to-surface-mapping "$cov_vol" "$midthick" "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii" \
      -ribbon-constrained "$white" "$pial" -volume-roi "$goodvoxels_mask"
    wb_command -metric-dilate "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii" "$midthick" 10 "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii" -nearest
    wb_command -metric-mask "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii" "$roi" "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii"
    wb_command -metric-resample "$gv_dir/$hemisphere.cov_after_mask.native.shape.gii" "$reg_native" "$reg_fsLR32k" \
      ADAP_BARY_AREA "$gv_dir/$hemisphere.cov_after_mask.32k_fs_LR.shape.gii" -area-surfs "$midthick" "$midthick_fsLR32k" -current-roi "$roi"
    wb_command -metric-dilate "$gv_dir/$hemisphere.cov_after_mask.32k_fs_LR.shape.gii" "$midthick_fsLR32k" 30 "$gv_dir/$hemisphere.cov_after_mask.32k_fs_LR.shape.gii" -nearest
    wb_command -metric-mask "$gv_dir/$hemisphere.cov_after_mask.32k_fs_LR.shape.gii" "$roi_fsLR32k" "$gv_dir/$hemisphere.cov_after_mask.32k_fs_LR.shape.gii"
  done

  wb_command -cifti-create-dense-scalar "$gv_dir/cov_before_mask.dscalar.nii" \
    -volume "$cov_vol" "$subcortical_rois" \
    -left-metric "$gv_dir/lh.cov_before_mask.32k_fs_LR.shape.gii" -roi-left "$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.L.atlasroi.32k_fs_LR.shape.gii" \
    -right-metric "$gv_dir/rh.cov_before_mask.32k_fs_LR.shape.gii" -roi-right "$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.R.atlasroi.32k_fs_LR.shape.gii"
  wb_command -set-map-names "$gv_dir/cov_before_mask.dscalar.nii" -map 1 "${subject}_cov_before_mask"

  wb_command -cifti-create-dense-scalar "$gv_dir/cov_after_mask.dscalar.nii" \
    -volume "$cov_vol" "$subcortical_sampling_mask" \
    -left-metric "$gv_dir/lh.cov_after_mask.32k_fs_LR.shape.gii" -roi-left "$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.L.atlasroi.32k_fs_LR.shape.gii" \
    -right-metric "$gv_dir/rh.cov_after_mask.32k_fs_LR.shape.gii" -roi-right "$subdir/anat/MNINonLinear/fsaverage_LR32k/$subject.R.atlasroi.32k_fs_LR.shape.gii"
  wb_command -set-map-names "$gv_dir/cov_after_mask.dscalar.nii" -map 1 "${subject}_cov_after_mask"
}

ensure_goodvoxels_sampling_mask() {
  local out_dir="$1"
  local ref_vol="$2"
  local cortical_ribbon_mask="$3"
  local subcortical_rois="$4"
  local subdir="$5"
  local subject="$6"
  local gv_dir="$out_dir/GoodVoxels"
  local goodvoxels_mask="$gv_dir/GoodVoxels.nii.gz"
  local cortical_sampling_mask="$gv_dir/CorticalRibbonGoodVoxelsMask.nii.gz"
  local subcortical_sampling_mask="$gv_dir/SubcorticalGoodVoxelsROIs.nii.gz"
  local subcortical_sampling_raw="$gv_dir/SubcorticalGoodVoxelsROIs_raw.nii.gz"
  local subcortical_label_table="$gv_dir/SubcorticalLabelTable.txt"
  local ribbon_only="$gv_dir/ribbon_only.nii.gz"

  mkdir -p "$gv_dir"
  rm -f "$gv_dir"/*

  cp "$cortical_ribbon_mask" "$ribbon_only"

  fslmaths "$ref_vol" -Tmean "$gv_dir/mean.nii.gz" -odt float
  fslmaths "$ref_vol" -Tstd "$gv_dir/std.nii.gz" -odt float
  fslmaths "$gv_dir/std.nii.gz" -div "$gv_dir/mean.nii.gz" "$gv_dir/cov.nii.gz"
  fslmaths "$gv_dir/cov.nii.gz" -mas "$ribbon_only" "$gv_dir/cov_ribbon.nii.gz"

  local cov_ribbon_mean cov_ribbon_mean_inv std mean upper
  cov_ribbon_mean="$(fslstats "$gv_dir/cov_ribbon.nii.gz" -M)"
  [[ -n "$cov_ribbon_mean" ]] || die "failed to compute GoodVoxels ribbon mean for $ref_vol"
  [[ "$cov_ribbon_mean" != "0" ]] || die "GoodVoxels ribbon mean is zero for $ref_vol"
  cov_ribbon_mean_inv="$(awk -v x="$cov_ribbon_mean" 'BEGIN { printf "%.10f", 1 / x }')"

  fslmaths "$gv_dir/cov_ribbon.nii.gz" -mul "$cov_ribbon_mean_inv" "$gv_dir/cov_ribbon_norm.nii.gz"
  fslmaths "$gv_dir/cov_ribbon_norm.nii.gz" -bin -s "$VOL2SURF_GOOD_VOXELS_SIGMA_MM" "$gv_dir/SmoothNorm.nii.gz"
  fslmaths "$gv_dir/cov_ribbon_norm.nii.gz" -s "$VOL2SURF_GOOD_VOXELS_SIGMA_MM" \
    -div "$gv_dir/SmoothNorm.nii.gz" -dilD "$gv_dir/cov_ribbon_norm_s${VOL2SURF_GOOD_VOXELS_SIGMA_MM}.nii.gz"
  fslmaths "$gv_dir/cov.nii.gz" -mul "$cov_ribbon_mean_inv" \
    -div "$gv_dir/cov_ribbon_norm_s${VOL2SURF_GOOD_VOXELS_SIGMA_MM}.nii.gz" "$gv_dir/cov_norm_modulate.nii.gz"
  fslmaths "$gv_dir/cov_norm_modulate.nii.gz" -mas "$ribbon_only" "$gv_dir/cov_norm_modulate_ribbon.nii.gz"

  std="$(fslstats "$gv_dir/cov_norm_modulate_ribbon.nii.gz" -S)"
  mean="$(fslstats "$gv_dir/cov_norm_modulate_ribbon.nii.gz" -M)"
  upper="$(awk -v mean="$mean" -v std="$std" -v factor="$VOL2SURF_GOOD_VOXELS_FACTOR" 'BEGIN { printf "%.10f", mean + (std * factor) }')"
  printf '%s\n' "$std" > "$gv_dir/std_scalar.txt"
  printf '%s\n' "$mean" > "$gv_dir/mean_scalar.txt"
  printf '%s\n' "$upper" > "$gv_dir/upper_scalar.txt"

  fslmaths "$gv_dir/mean.nii.gz" -bin "$gv_dir/mask.nii.gz"
  fslmaths "$gv_dir/cov_norm_modulate.nii.gz" -thr "$upper" -bin -sub "$gv_dir/mask.nii.gz" -mul -1 "$goodvoxels_mask"

  fslmaths "$goodvoxels_mask" -mas "$ribbon_only" "$cortical_sampling_mask"
  wb_command -volume-label-export-table "$subcortical_rois" 1 "$subcortical_label_table"
  fslmaths "$subcortical_rois" -mas "$goodvoxels_mask" "$subcortical_sampling_raw"
  wb_command -volume-label-import "$subcortical_sampling_raw" "$subcortical_label_table" "$subcortical_sampling_mask" \
    -discard-others -unlabeled-value 0 -drop-unused-labels
  write_goodvoxels_surface_qc "$gv_dir" "$gv_dir/cov.nii.gz" "$goodvoxels_mask" "$subcortical_rois" "$subcortical_sampling_mask" "$subdir" "$subject"
  cleanup_goodvoxels_workspace "$gv_dir" "$goodvoxels_mask" "$cortical_sampling_mask" "$subcortical_sampling_mask"
  printf '%s\n%s\n' "$cortical_sampling_mask" "$subcortical_sampling_mask"
}

# Parse comma-separated VOL2SURF_INPUTS.
mapfile -t CIFTI_TAGS < <(echo "$InputSpec" | tr ',' '\n' | sed 's/^ *//; s/ *$//' | awk 'NF')
[[ "${#CIFTI_TAGS[@]}" -gt 0 ]] || die "no valid tags parsed from InputSpec: '$InputSpec'"

mapfile -t SESSION_DIRS < <(
  find -L "$Subdir/func/$FuncDirName" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V
)

for SES_DIR in "${SESSION_DIRS[@]}"; do
  s="${SES_DIR##*/}"
  s="${s#session_}"
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  (( s >= StartSession )) || continue

  mapfile -t RUN_DIRS < <(
    find -L "$SES_DIR" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V
  )

  for OUT_DIR in "${RUN_DIRS[@]}"; do
    r="${OUT_DIR##*/}"
    r="${r#run_}"
    [[ "$r" =~ ^[0-9]+$ ]] || continue

    RunGoodVoxelsReady=0
    RunVolumeRoiMask=""
    RunCiftiVolumeMask=""
    if [[ "$VOL2SURF_USE_GOOD_VOXELS_MASK" == "1" ]]; then
      [[ "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" == "1" ]] || die "GoodVoxels masking requires VOL2SURF_USE_CORTICAL_RIBBON_MASK=1"
      RunGoodVoxelsRefNii="$OUT_DIR/${FuncFilePrefix}_OCME.nii.gz"
      [[ -f "$RunGoodVoxelsRefNii" ]] || die "missing GoodVoxels reference volume: $RunGoodVoxelsRefNii"
    fi

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
      subcortical_rois="$Subdir/func/rois/Subcortical_ROIs_acpc.nii.gz"

      CorticalRibbonMask=""
      if [[ "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" == "1" ]]; then
        CorticalRibbonMask="$(ensure_cortical_ribbon_func_mask "$Subdir" "$in_nii" "$MEDIR" "$AtlasSpace")"
      fi

      VolumeRoiMask=""
      CiftiVolumeMask="$subcortical_rois"
      CiftiVolumeData="$in_nii"
      if [[ "$VOL2SURF_USE_CORTICAL_RIBBON_MASK" == "1" ]]; then
        VolumeRoiMask="$CorticalRibbonMask"
      fi
      if [[ "$VOL2SURF_USE_GOOD_VOXELS_MASK" == "1" ]]; then
        if [[ "$RunGoodVoxelsReady" == "0" ]]; then
          [[ -n "$CorticalRibbonMask" ]] || die "GoodVoxels masking requires VOL2SURF_USE_CORTICAL_RIBBON_MASK=1"
          mapfile -t goodvoxels_masks < <(ensure_goodvoxels_sampling_mask "$OUT_DIR" "$RunGoodVoxelsRefNii" "$CorticalRibbonMask" "$subcortical_rois" "$Subdir" "$Subject")
          [[ "${#goodvoxels_masks[@]}" -eq 2 ]] || die "expected cortical and subcortical GoodVoxels masks for $RunGoodVoxelsRefNii"
          RunVolumeRoiMask="${goodvoxels_masks[0]}"
          RunCiftiVolumeMask="${goodvoxels_masks[1]}"
          RunGoodVoxelsReady=1
        fi
        VolumeRoiMask="$RunVolumeRoiMask"
        # Keep a fixed subcortical ROI across runs so downstream concat sees a stable grayordinate axis.
        # Apply the run-specific GoodVoxels mask to the volume data instead of shrinking the CIFTI axis.
        CiftiVolumeMask="$subcortical_rois"
      fi
      UseGoodVoxelsFill="$VOL2SURF_USE_GOOD_VOXELS_MASK"
      if [[ "$UseGoodVoxelsFill" == "1" ]]; then
        CiftiVolumeData="$OUT_DIR/${FuncFilePrefix}_${tag}${out_suffix}_subcort_goodvox.nii.gz"
        fslmaths "$in_nii" -mas "$RunCiftiVolumeMask" "$CiftiVolumeData"
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

        if [[ -n "$VolumeRoiMask" ]]; then
          wb_command -volume-to-surface-mapping "$in_nii" "$MIDTHICK" "$OUT_DIR/$hemisphere.native.shape.gii" \
            -ribbon-constrained "$WHITE" "$PIAL" -volume-roi "$VolumeRoiMask"
        else
          wb_command -volume-to-surface-mapping "$in_nii" "$MIDTHICK" "$OUT_DIR/$hemisphere.native.shape.gii" \
            -ribbon-constrained "$WHITE" "$PIAL"
        fi

        wb_command -metric-dilate "$OUT_DIR/$hemisphere.native.shape.gii" "$MIDTHICK" 10 "$OUT_DIR/$hemisphere.native.shape.gii" -nearest
        wb_command -metric-mask "$OUT_DIR/$hemisphere.native.shape.gii" "$ROI" "$OUT_DIR/$hemisphere.native.shape.gii"
        wb_command -metric-resample "$OUT_DIR/$hemisphere.native.shape.gii" "$REG_MSMSulc" "$REG_MSMSulc_FSLR32k" \
          ADAP_BARY_AREA "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" -area-surfs "$MIDTHICK" "$MIDTHICK_FSLR32k" -current-roi "$ROI"
        if [[ "$UseGoodVoxelsFill" == "1" ]]; then
          wb_command -metric-dilate "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" "$MIDTHICK_FSLR32k" 30 "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" -nearest
        fi
        wb_command -metric-mask "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii" "$ROI_FSLR32k" "$OUT_DIR/$hemisphere.32k_fs_LR.shape.gii"
      done

      tr=$(cat "$OUT_DIR/TR.txt")
      if [[ "$UseGoodVoxelsFill" == "1" ]]; then
        out_dt_tmp="$OUT_DIR/${FuncFilePrefix}_${tag}${out_suffix}_raw.dtseries.nii"
        wb_command -cifti-create-dense-timeseries "$out_dt_tmp" -volume "$CiftiVolumeData" "$CiftiVolumeMask" \
          -left-metric "$OUT_DIR/lh.32k_fs_LR.shape.gii" -roi-left "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.L.atlasroi.32k_fs_LR.shape.gii" \
          -right-metric "$OUT_DIR/rh.32k_fs_LR.shape.gii" -roi-right "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.R.atlasroi.32k_fs_LR.shape.gii" \
          -timestep "$tr"
        wb_command -cifti-dilate "$out_dt_tmp" COLUMN 0 30 "$out_dt" \
          -left-surface "$Subdir/anat/T1w/fsaverage_LR32k/$Subject.L.midthickness.32k_fs_LR.surf.gii" \
          -right-surface "$Subdir/anat/T1w/fsaverage_LR32k/$Subject.R.midthickness.32k_fs_LR.surf.gii"
        rm -f "$out_dt_tmp" "$CiftiVolumeData"
      else
        wb_command -cifti-create-dense-timeseries "$out_dt" -volume "$CiftiVolumeData" "$CiftiVolumeMask" \
          -left-metric "$OUT_DIR/lh.32k_fs_LR.shape.gii" -roi-left "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.L.atlasroi.32k_fs_LR.shape.gii" \
          -right-metric "$OUT_DIR/rh.32k_fs_LR.shape.gii" -roi-right "$Subdir/anat/MNINonLinear/fsaverage_LR32k/$Subject.R.atlasroi.32k_fs_LR.shape.gii" \
          -timestep "$tr"
      fi

      rm -f "$OUT_DIR"/*shape*
    done
  done
done

#!/usr/bin/env bash
# MEDIC-aware functional-to-anatomical coregistration.
#
# This module assumes per-run MEDIC displacement maps already exist under:
#   func/<task>/session_<s>/run_<r>/MEDIC/
# It builds MEDIC-unwarped SBRefs, estimates the unwarped SBRef -> T1w_acpc
# registration, and writes the run-level pointer files consumed by
# mefmri_func_headmotion_medic.sh.

set -euo pipefail
shopt -s nullglob

if [[ "$#" -lt 7 || "$#" -gt 10 ]]; then
  echo "Usage: mefmri_func_coreg_medic.sh <MEDIR> <Subject> <StudyFolder> <AtlasTemplate> <DOF> <NTHREADS> <StartSession> [AtlasSpace] [FuncDirName] [FuncFilePrefix]" >&2
  exit 2
fi

MEDIR="$1"
Subject="$2"
StudyFolder="$3"
AtlasTemplate="$4"
DOF="$5"
NTHREADS="$6"
StartSession="$7"
AtlasSpace="${8:-${AtlasSpace:-T1w}}"
FuncDirName="${9:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${10:-${FUNC_FILE_PREFIX:-Rest}}"
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-$FuncDirName}"

[[ "${MEDIR:0:1}" == "/" ]] || MEDIR="$(cd "$MEDIR" && pwd)"
[[ "${StudyFolder:0:1}" == "/" ]] || StudyFolder="$(cd "$StudyFolder" && pwd)"
[[ "${AtlasTemplate:0:1}" == "/" ]] || AtlasTemplate="$(cd "$(dirname "$AtlasTemplate")" && pwd)/$(basename "$AtlasTemplate")"

Subdir="$StudyFolder/$Subject"
FuncRoot="$Subdir/func/$FuncDirName"
UnprocRoot="$Subdir/func/unprocessed/$FuncDirName"
XfmsDir="$Subdir/func/xfms/$FuncXfmsDir"
QC_OUTPUT_DIRNAME="${QC_OUTPUT_DIRNAME:-QC}"
COREG_QA_DIR="${COREG_QA_DIR:-$FuncRoot/${QC_OUTPUT_DIRNAME}/Coreg}"
COREG_PYTHON="${COREG_PYTHON:-${PIPELINE_PYTHON:-python3}}"
MedicReferenceModule="${FUNC_MEDIC_REFERENCE_MODULE:-$MEDIR/modules/mefmri_func_medic_reference.sh}"
MEDIC_REFERENCE_POLICY="${MEDIC_REFERENCE_POLICY:-first}"

case "$AtlasSpace" in
  T1w|MNINonlinear) ;;
  *)
    echo "ERROR: mefmri_func_coreg_medic.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)" >&2
    exit 2
    ;;
esac

die() {
  echo "ERROR: $*" >&2
  exit 2
}

MedicPreserveDir=""

preserve_run_local_medic_dirs() {
  MedicPreserveDir="$(mktemp -d -p /tmp "mefmri_medic_preserve_XXXXXX")"
  [[ -d "$FuncRoot" ]] || return 0

  local medic_dir run_dir session_dir rel_dir dest_dir
  while IFS= read -r medic_dir; do
    run_dir="$(dirname "$medic_dir")"
    session_dir="$(dirname "$run_dir")"
    rel_dir="$(basename "$session_dir")/$(basename "$run_dir")"
    dest_dir="$MedicPreserveDir/$rel_dir"
    mkdir -p "$dest_dir"
    mv "$medic_dir" "$dest_dir/MEDIC"
  done < <(find -L "$FuncRoot" -mindepth 3 -maxdepth 3 -type d -name MEDIC | sort -V)
}

restore_run_local_medic_dirs() {
  [[ -n "$MedicPreserveDir" && -d "$MedicPreserveDir" ]] || return 0

  local preserved rel_dir dest_dir
  while IFS= read -r preserved; do
    rel_dir="${preserved#"$MedicPreserveDir"/}"
    rel_dir="${rel_dir%/MEDIC}"
    dest_dir="$FuncRoot/$rel_dir"
    mkdir -p "$dest_dir"
    rm -rf "$dest_dir/MEDIC"
    mv "$preserved" "$dest_dir/MEDIC"
  done < <(find -L "$MedicPreserveDir" -mindepth 3 -maxdepth 3 -type d -name MEDIC | sort -V)

  rm -rf "$MedicPreserveDir"
  MedicPreserveDir=""
}

ensure_white_wmseg() {
  local fs_mri_dir="$Subdir/anat/T1w/$Subject/mri"
  local aseg_mgz="$fs_mri_dir/aparc+aseg.mgz"
  local white_mgz="$fs_mri_dir/white.mgz"
  local white_nii="$fs_mri_dir/white.nii.gz"
  if [[ -f "$white_nii" ]]; then
    echo "$white_nii"
    return 0
  fi
  [[ -f "$aseg_mgz" ]] || die "missing FreeSurfer aseg for WM segmentation: $aseg_mgz"
  mri_binarize --i "$aseg_mgz" --wm --o "$white_mgz" >/dev/null 2>&1
  mri_convert -i "$white_mgz" -o "$white_nii" --like "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz" >/dev/null 2>&1
  [[ -f "$white_nii" ]] || die "failed to create WM segmentation: $white_nii"
  echo "$white_nii"
}

ensure_white_deformed() {
  export SUBJECTS_DIR="$Subdir/anat/T1w"
  local fs_subj="freesurfer"
  local t1acpc="$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz"
  local origmgz="$SUBJECTS_DIR/$fs_subj/mri/orig.mgz"
  local lh_def="$SUBJECTS_DIR/$fs_subj/surf/lh.white.deformed"
  local rh_def="$SUBJECTS_DIR/$fs_subj/surf/rh.white.deformed"

  rm -rf "$SUBJECTS_DIR/$fs_subj"
  cp -rf "$Subdir/anat/T1w/$Subject" "$SUBJECTS_DIR/$fs_subj"

  if [[ -f "$lh_def" && -f "$rh_def" ]]; then
    echo "[medic-coreg] Found existing white.deformed surfaces"
    return 0
  fi

  [[ -f "$origmgz" ]] || die "missing $origmgz"
  [[ -f "$t1acpc" ]] || die "missing $t1acpc"
  [[ -f "$SUBJECTS_DIR/$fs_subj/surf/lh.white" ]] || die "missing lh.white in $SUBJECTS_DIR/$fs_subj/surf"
  [[ -f "$SUBJECTS_DIR/$fs_subj/surf/rh.white" ]] || die "missing rh.white in $SUBJECTS_DIR/$fs_subj/surf"

  local reg_tmp
  reg_tmp="$(mktemp -p /tmp "${fs_subj}_orig2acpc_XXXXXX.dat")"
  tkregister2 --mov "$t1acpc" --targ "$origmgz" --noedit --regheader --reg "$reg_tmp"
  mri_surf2surf --s "$fs_subj" --hemi lh --sval-xyz white --reg "$reg_tmp" --tval-xyz "$t1acpc" --tval "$lh_def"
  mri_surf2surf --s "$fs_subj" --hemi rh --sval-xyz white --reg "$reg_tmp" --tval-xyz "$t1acpc" --tval "$rh_def"
  rm -f "$reg_tmp"

  [[ -f "$lh_def" && -f "$rh_def" ]] || die "failed to create white.deformed surfaces"
}

write_grid_anatomicals_and_masks() {
  mkdir -p "$XfmsDir"

  if [[ -f "$Subdir/anat/T1w/T2w_acpc_dc_restore.nii.gz" ]]; then
    flirt -interp nearestneighbour -in "$Subdir/anat/T1w/T2w_acpc_dc_restore.nii.gz" -ref "$AtlasTemplate" -out "$XfmsDir/T2w_acpc_func.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
  fi
  if [[ -f "$Subdir/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz" ]]; then
    flirt -interp nearestneighbour -in "$Subdir/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz" -ref "$AtlasTemplate" -out "$XfmsDir/T2w_acpc_brain_func.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
  fi

  flirt -interp nearestneighbour -in "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz" -ref "$AtlasTemplate" -out "$XfmsDir/T1w_acpc_func.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
  flirt -interp nearestneighbour -in "$Subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" -ref "$AtlasTemplate" -out "$XfmsDir/T1w_acpc_brain_func.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
  fslmaths "$XfmsDir/T1w_acpc_brain_func.nii.gz" -bin "$XfmsDir/T1w_acpc_brain_func_mask.nii.gz"

  local ribbon_src="$Subdir/anat/T1w/CorticalRibbon.nii.gz"
  if [[ ! -f "$ribbon_src" && -f "$Subdir/anat/T1w/CorticalRibbon.ni.gz" ]]; then
    ribbon_src="$Subdir/anat/T1w/CorticalRibbon.ni.gz"
  fi
  if [[ ! -f "$ribbon_src" ]]; then
    local lh_ribbon="$Subdir/anat/T1w/$Subject/mri/lh.ribbon.mgz"
    local rh_ribbon="$Subdir/anat/T1w/$Subject/mri/rh.ribbon.mgz"
    if [[ -f "$lh_ribbon" && -f "$rh_ribbon" ]]; then
      local tmp_lh="$XfmsDir/.tmp_lh.ribbon.nii.gz"
      local tmp_rh="$XfmsDir/.tmp_rh.ribbon.nii.gz"
      ribbon_src="$XfmsDir/.tmp_CorticalRibbon_auto.nii.gz"
      mri_convert -i "$lh_ribbon" -o "$tmp_lh" --like "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz"
      mri_convert -i "$rh_ribbon" -o "$tmp_rh" --like "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz"
      fslmaths "$tmp_lh" -add "$tmp_rh" "$ribbon_src"
      fslmaths "$ribbon_src" -bin "$ribbon_src"
      rm -f "$tmp_lh" "$tmp_rh"
    else
      echo "[medic-coreg] WARNING: missing cortical ribbon; falling back to T1w brain mask"
      ribbon_src="$Subdir/anat/T1w/T1w_acpc_brain_mask.nii.gz"
    fi
  fi
  flirt -interp nearestneighbour -in "$ribbon_src" -ref "$AtlasTemplate" -out "$XfmsDir/CorticalRibbon_acpc_func_mask.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
  fslmaths "$XfmsDir/CorticalRibbon_acpc_func_mask.nii.gz" -bin "$XfmsDir/CorticalRibbon_acpc_func_mask.nii.gz"

  local nonlin_warp="$Subdir/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz"
  if [[ -f "$nonlin_warp" ]]; then
    applywarp --interp=nn --in="$ribbon_src" --ref="$AtlasTemplate" --warp="$nonlin_warp" --out="$XfmsDir/CorticalRibbon_nonlin_func_mask.nii.gz"
    fslmaths "$XfmsDir/CorticalRibbon_nonlin_func_mask.nii.gz" -bin "$XfmsDir/CorticalRibbon_nonlin_func_mask.nii.gz"
    [[ -f "$Subdir/anat/MNINonLinear/T1w_restore_brain.nii.gz" ]] && flirt -interp nearestneighbour -in "$Subdir/anat/MNINonLinear/T1w_restore_brain.nii.gz" -ref "$AtlasTemplate" -out "$XfmsDir/T1w_nonlin_brain_func.nii.gz" -applyxfm -init "$MEDIR/res0urces/ident.mat"
    [[ -f "$XfmsDir/T1w_nonlin_brain_func.nii.gz" ]] && fslmaths "$XfmsDir/T1w_nonlin_brain_func.nii.gz" -bin "$XfmsDir/T1w_nonlin_brain_func_mask.nii.gz"
  elif [[ "$AtlasSpace" == "MNINonlinear" ]]; then
    die "missing nonlinear warp required for AtlasSpace=MNINonlinear: $nonlin_warp"
  fi
  rm -f "$XfmsDir/.tmp_CorticalRibbon_auto.nii.gz"
}

[[ -d "$UnprocRoot" ]] || die "missing unprocessed functional root: $UnprocRoot"
[[ -x "$MedicReferenceModule" || -f "$MedicReferenceModule" ]] || die "missing MEDIC reference module: $MedicReferenceModule"

echo "[medic-coreg] Subject=$Subject FuncDirName=$FuncDirName AtlasSpace=$AtlasSpace"
echo "[medic-coreg] Building MEDIC-unwarped SBRefs using policy=$MEDIC_REFERENCE_POLICY"

preserve_run_local_medic_dirs
trap restore_run_local_medic_dirs EXIT
"$COREG_PYTHON" "$MEDIR/lib/find_epi_params.py" \
  --subdir "$Subdir" --func-name "$FuncDirName" --func-prefix "$FuncFilePrefix" --start-session "$StartSession" --qc-dirname "$QC_OUTPUT_DIRNAME" --no-fieldmap-mode
restore_run_local_medic_dirs
trap - EXIT

bash "$MedicReferenceModule" "$MEDIR" "$Subject" "$StudyFolder" "$NTHREADS" "$StartSession" "$FuncDirName" "$FuncFilePrefix" "$MEDIC_REFERENCE_POLICY"

WMSEG_NII="$(ensure_white_wmseg)"
mkdir -p "$XfmsDir" "$COREG_QA_DIR" "$FuncRoot/AverageSBref"
WDIR="$FuncRoot/AverageSBref"
rm -f "$WDIR"/MEDICSBref_*.nii.gz "$WDIR"/TMP_MEDIC*.nii.gz

mapfile -t SessionDirs < <(find -L "$FuncRoot" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V)
for SessionDir in "${SessionDirs[@]}"; do
  s="${SessionDir##*/session_}"
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  (( s >= StartSession )) || continue
  mapfile -t RunDirs < <(find -L "$SessionDir" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
  for RunDir in "${RunDirs[@]}"; do
    r="${RunDir##*/run_}"
    [[ "$r" =~ ^[0-9]+$ ]] || continue
    ref="$RunDir/MEDIC/reference/SBref_unwarped.nii.gz"
    [[ -f "$ref" ]] || die "missing MEDIC-unwarped SBRef: $ref"
    cp "$ref" "$WDIR/MEDICSBref_${s}_${r}.nii.gz"
  done
done

images=("$WDIR"/MEDICSBref_*.nii.gz)
[[ "${#images[@]}" -gt 0 ]] || die "no MEDIC-unwarped SBRefs were found"
if [[ "${#images[@]}" -gt 1 ]]; then
  "$MEDIR/res0urces/FuncAverage" -n -o "$XfmsDir/AvgMEDICSBref.nii.gz" "${images[@]}"
else
  cp "${images[0]}" "$XfmsDir/AvgMEDICSBref.nii.gz"
fi

ensure_white_deformed

"$MEDIR/res0urces/epi_reg_dof" --dof="$DOF" \
  --epi="$XfmsDir/AvgMEDICSBref.nii.gz" \
  --t1="$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz" \
  --t1brain="$Subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
  --out="$XfmsDir/AvgMEDICSBref2acpc_EpiReg" \
  --wmseg="$WMSEG_NII"

convertwarp --ref="$AtlasTemplate" --premat="$XfmsDir/AvgMEDICSBref2acpc_EpiReg.mat" --out="$XfmsDir/AvgMEDICSBref2acpc_EpiReg_warp.nii.gz"
applywarp --interp=spline --in="$XfmsDir/AvgMEDICSBref.nii.gz" --ref="$AtlasTemplate" --out="$XfmsDir/AvgMEDICSBref2acpc_EpiReg.nii.gz" --warp="$XfmsDir/AvgMEDICSBref2acpc_EpiReg_warp.nii.gz"

bbregister --s freesurfer --mov "$XfmsDir/AvgMEDICSBref2acpc_EpiReg.nii.gz" --init-reg "$MEDIR/res0urces/eye.dat" --surf white.deformed --bold --reg "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.dat" --6 --o "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.nii.gz"
tkregister2 --s freesurfer --noedit --reg "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.dat" --mov "$XfmsDir/AvgMEDICSBref2acpc_EpiReg.nii.gz" --targ "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz" --fslregout "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.mat"

convertwarp --warp1="$XfmsDir/AvgMEDICSBref2acpc_EpiReg_warp.nii.gz" --postmat="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.mat" --ref="$AtlasTemplate" --out="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_warp.nii.gz"
applywarp --interp=spline --in="$XfmsDir/AvgMEDICSBref.nii.gz" --ref="$AtlasTemplate" --out="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.nii.gz" --warp="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_warp.nii.gz"
invwarp --ref="$XfmsDir/AvgMEDICSBref.nii.gz" -w "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_warp.nii.gz" -o "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_inv_warp.nii.gz"
convert_xfm -omat "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_full.mat" -concat "$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR.mat" "$XfmsDir/AvgMEDICSBref2acpc_EpiReg.mat"

if [[ -f "$Subdir/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz" ]]; then
  convertwarp --ref="$AtlasTemplate" --warp1="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_warp.nii.gz" --warp2="$Subdir/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz" --out="$XfmsDir/AvgMEDICSBref2nonlin_EpiReg+BBR_warp.nii.gz"
  applywarp --interp=spline --in="$XfmsDir/AvgMEDICSBref.nii.gz" --ref="$AtlasTemplate" --out="$XfmsDir/AvgMEDICSBref2nonlin_EpiReg+BBR.nii.gz" --warp="$XfmsDir/AvgMEDICSBref2nonlin_EpiReg+BBR_warp.nii.gz"
fi

write_grid_anatomicals_and_masks

if [[ "$AtlasSpace" == "MNINonlinear" ]]; then
  PointerWarp="$XfmsDir/AvgMEDICSBref2nonlin_EpiReg+BBR_warp.nii.gz"
else
  PointerWarp="$XfmsDir/AvgMEDICSBref2acpc_EpiReg+BBR_warp.nii.gz"
fi
[[ -f "$PointerWarp" ]] || die "missing selected MEDIC coreg warp: $PointerWarp"

PointerLog="$COREG_QA_DIR/MEDICCoregPointerSelection.tsv"
echo -e "session\trun\ttarget\twarp\tdisplacement\tphase_encoding_axis" > "$PointerLog"

for SessionDir in "${SessionDirs[@]}"; do
  s="${SessionDir##*/session_}"
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  (( s >= StartSession )) || continue
  mapfile -t RunDirs < <(find -L "$SessionDir" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
  for RunDir in "${RunDirs[@]}"; do
    r="${RunDir##*/run_}"
    [[ "$r" =~ ^[0-9]+$ ]] || continue
    dmap="$RunDir/MEDIC/${FuncFilePrefix}_S${s}_R${r}_displacementmaps.nii"
    ref="$RunDir/MEDIC/reference/SBref_unwarped.nii.gz"
    policy="$RunDir/MEDIC/reference/reference_policy.tsv"
    ped="$(awk -F '\t' '$1=="phase_encoding_direction"{print $2}' "$policy" 2>/dev/null)"
    [[ -f "$dmap" ]] || die "missing MEDIC displacement map: $dmap"
    [[ -f "$ref" ]] || die "missing MEDIC-unwarped SBRef: $ref"
    [[ -n "$ped" ]] || die "could not read phase_encoding_direction from $policy"

    echo "$XfmsDir/AvgMEDICSBref.nii.gz" > "$RunDir/IntermediateCoregTarget.txt"
    echo "$PointerWarp" > "$RunDir/Intermediate2ACPCWarp.txt"
    rm -f "$RunDir/Intermediate2ACPCMat.txt"
    echo "$dmap" > "$RunDir/MEDICDisplacementMap.txt"
    echo "$ped" > "$RunDir/MEDICPhaseEncodingAxis.txt"
    echo "$ref" > "$RunDir/MEDICUnwarpedSBRef.txt"
    echo -e "${s}\t${r}\t$XfmsDir/AvgMEDICSBref.nii.gz\t$PointerWarp\t$dmap\t$ped" >> "$PointerLog"
  done
done

"$MEDIR/lib/make_precise_subcortical_labels.sh" "$Subdir" "$AtlasTemplate" "$MEDIR"
rm -rf "$Subdir/anat/T1w/freesurfer"

echo "[medic-coreg] complete"

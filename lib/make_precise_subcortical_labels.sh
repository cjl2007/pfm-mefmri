#!/usr/bin/env bash
# No-MATLAB replacement for make_precise_subcortical_labels.m
# Usage: make_precise_subcortical_labels.sh <Subdir> <AtlasTemplate> <MEDIR>

set -euo pipefail
IFS=$'\n\t'

Subdir="${1:?missing Subdir}"
AtlasTemplate="${2:?missing AtlasTemplate}"
MEDIR="${3:?missing MEDIR}"

LABELS=(8 47 26 58 18 54 11 50 17 53 13 52 12 51 10 49 16 28 60)
IDENT_MAT="${MEDIR}/res0urces/ident.mat"
SUBCORT_LABELS_TXT="${SUBCORT_LABELS_TXT:-${SUBCORTICAL_LABELS_TXT:-${MEDIR}/res0urces/SubcorticalLabels.txt}}"
ROIS_DIR="${Subdir}/func/rois"
TMP_DIR="${ROIS_DIR}/tmp"

mkdir -p "${ROIS_DIR}" "${TMP_DIR}"

need_file() { [[ -f "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 2; }; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 2; }; }

need_cmd mri_binarize
need_cmd flirt
need_cmd fslmaths
need_cmd fslmerge
need_cmd wb_command

need_file "${IDENT_MAT}"
need_file "${SUBCORT_LABELS_TXT}"
need_file "${AtlasTemplate}"

make_space_labels() {
  local aparc="$1"
  local out="$2"

  need_file "${aparc}"
  rm -rf "${TMP_DIR}"
  mkdir -p "${TMP_DIR}"

  local idx=0
  for label in "${LABELS[@]}"; do
    idx=$((idx + 1))
    mri_binarize --i "${aparc}" --match "${label}" --o "${TMP_DIR}/Label${idx}.nii.gz" >/dev/null 2>&1
    flirt -interp nearestneighbour -in "${TMP_DIR}/Label${idx}.nii.gz" -ref "${AtlasTemplate}" \
      -applyxfm -init "${IDENT_MAT}" -out "${TMP_DIR}/Label${idx}_Interp.nii.gz" >/dev/null 2>&1
    fslmaths "${TMP_DIR}/Label${idx}_Interp.nii.gz" -mul "${label}" "${TMP_DIR}/Label${idx}_Final.nii.gz" >/dev/null 2>&1
  done

  fslmerge -t "${TMP_DIR}/FinalLabels.nii.gz" "${TMP_DIR}"/Label*_Final.nii.gz >/dev/null 2>&1
  fslmaths "${TMP_DIR}/FinalLabels.nii.gz" -Tmax "${TMP_DIR}/FinalLabels.nii.gz" >/dev/null 2>&1
  wb_command -volume-label-import "${TMP_DIR}/FinalLabels.nii.gz" "${SUBCORT_LABELS_TXT}" "${out}" -discard-others >/dev/null 2>&1
}

make_space_labels "${Subdir}/anat/T1w/aparc+aseg.nii.gz" "${ROIS_DIR}/Subcortical_ROIs_acpc.nii.gz"
make_space_labels "${Subdir}/anat/MNINonLinear/aparc+aseg.nii.gz" "${ROIS_DIR}/Subcortical_ROIs_nonlin.nii.gz"

rm -rf "${TMP_DIR}"

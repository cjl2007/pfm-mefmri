#!/usr/bin/env bash
# Apply manual Infomap review-table edits without rerunning Infomap.
set -euo pipefail

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
MEDIR="${3:?missing MEDIR}"
_unused_start_session="${4:?missing StartSession}"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
_unused_prefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

SubjectDir="${StudyFolder}/${Subject}"
PFM_PYTHON="${PFM_PYTHON:-${PIPELINE_PYTHON:-python3}}"
PFM_ROOT_DIR="${PFM_OUTDIR:-${SubjectDir}/func/${FuncDirName}/PFM}"
PFM_INFOMAP_DIR="${PFM_ROOT_DIR}/Infomap"

PFM_PRIORS_MAT="${PFM_PRIORS_MAT:-${NETWORK_PRIORS_MAT:-${MEDIR}/res0urces/Priors.mat}}"
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="${PFM_INFOMAP_GRAPH_DENSITIES_EXPR:-0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001}"
PFM_INFOMAP_LABEL_DENSITY_INDEX="${PFM_INFOMAP_LABEL_DENSITY_INDEX:--1}"
PFM_INFOMAP_LABEL_UNASSIGNED_VALUE="${PFM_INFOMAP_LABEL_UNASSIGNED_VALUE:-21}"
PFM_INFOMAP_LABEL_WB_COMMAND="${PFM_INFOMAP_LABEL_WB_COMMAND:-wb_command}"
PFM_INFOMAP_UPDATE_TABLE_GLOB="${PFM_INFOMAP_UPDATE_TABLE_GLOB:-GraphDensity_*/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv}"
PFM_INFOMAP_UPDATE_OUTFILE="${PFM_INFOMAP_UPDATE_OUTFILE:-InfomapNetworkLabels_ManualAdjusted}"
PFM_INPUT_CIFTI="${PFM_INPUT_CIFTI:-}"
PFM_INPUT_TAG="${PFM_INPUT_TAG:-${CONCAT_INPUT_TAG:-OCME+MEICA+MGTR}}"
PFM_CONCAT_OUT_SUBDIR="${PFM_CONCAT_OUT_SUBDIR:-${CONCAT_OUT_SUBDIR:-ConcatenatedCiftis}}"
PFM_FD_THRESHOLD="${PFM_FD_THRESHOLD:-${CONCAT_FD_THRESHOLD:-0.3}}"
PFM_HOMOGENEITY_TEST_ENABLE="${PFM_HOMOGENEITY_TEST_ENABLE:-1}"
PFM_HOMOGENEITY_ROTATIONS_CIFTI="${PFM_HOMOGENEITY_ROTATIONS_CIFTI:-${MEDIR}/res0urces/Rotated_inds.dtseries.nii}"
PFM_HOMOGENEITY_N_ROTATIONS="${PFM_HOMOGENEITY_N_ROTATIONS:-1000}"
PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE="${PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE:-5}"
PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY="${PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY:-1000}"
PFM_HOMOGENEITY_ALPHA="${PFM_HOMOGENEITY_ALPHA:-0.05}"

L_MID="${SubjectDir}/anat/T1w/fsaverage_LR32k/${Subject}.L.midthickness.32k_fs_LR.surf.gii"
R_MID="${SubjectDir}/anat/T1w/fsaverage_LR32k/${Subject}.R.midthickness.32k_fs_LR.surf.gii"
COMMUNITIES="${PFM_INFOMAP_DIR}/Bipartite_PhysicalCommunities.dtseries.nii"
if [[ -z "$PFM_INPUT_CIFTI" ]]; then
  FD_TAG="${PFM_FD_THRESHOLD//./p}"
  PFM_INPUT_CIFTI="${SubjectDir}/func/${FuncDirName}/${PFM_CONCAT_OUT_SUBDIR}/${_unused_prefix}_${PFM_INPUT_TAG}_Concatenated+FDlt${FD_TAG}.dtseries.nii"
fi

[[ -d "$PFM_INFOMAP_DIR" ]] || { echo "ERROR: missing Infomap output dir: $PFM_INFOMAP_DIR"; exit 2; }
[[ -f "$COMMUNITIES" ]] || { echo "ERROR: missing Infomap communities CIFTI: $COMMUNITIES"; exit 2; }
[[ -f "$PFM_PRIORS_MAT" ]] || { echo "ERROR: missing Priors.mat: $PFM_PRIORS_MAT"; exit 2; }

echo "[pfm_update] scanning ${PFM_INFOMAP_DIR}/${PFM_INFOMAP_UPDATE_TABLE_GLOB}"
"$PFM_PYTHON" "$MEDIR/lib/pfm_infomap_manual_labels.py" \
  --communities-cifti "$COMMUNITIES" \
  --manual-corrections-glob "$PFM_INFOMAP_UPDATE_TABLE_GLOB" \
  --priors-mat "$PFM_PRIORS_MAT" \
  --outdir "$PFM_INFOMAP_DIR" \
  --outfile-prefix "$PFM_INFOMAP_UPDATE_OUTFILE" \
  --density-index "$PFM_INFOMAP_LABEL_DENSITY_INDEX" \
  --density-values "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR" \
  --unassigned-value "$PFM_INFOMAP_LABEL_UNASSIGNED_VALUE" \
  --left-surf "$L_MID" \
  --right-surf "$R_MID" \
  --wb-command "$PFM_INFOMAP_LABEL_WB_COMMAND" \
  --density-output-mode subdirs

ADJUSTED_CONSENSUS="${PFM_INFOMAP_DIR}/${PFM_INFOMAP_UPDATE_OUTFILE}_ModeConsensus.dlabel.nii"
if [[ "$PFM_HOMOGENEITY_TEST_ENABLE" == "1" && -f "$PFM_INPUT_CIFTI" && -f "$ADJUSTED_CONSENSUS" ]]; then
  HOMOG_ROT_ARGS=()
  if [[ -f "$PFM_HOMOGENEITY_ROTATIONS_CIFTI" ]]; then
    HOMOG_ROT_ARGS=( --rotations-cifti "$PFM_HOMOGENEITY_ROTATIONS_CIFTI" --n-rotations "$PFM_HOMOGENEITY_N_ROTATIONS" )
  else
    HOMOG_ROT_ARGS=( --n-rotations 0 )
  fi
  "$PFM_PYTHON" "$MEDIR/lib/pfm_homogeneity_test.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --labels-cifti "$ADJUSTED_CONSENSUS" \
    --outdir "${PFM_INFOMAP_DIR}/HomogeneityTest/ManualAdjustedConsensus" \
    --min-community-size "$PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE" \
    --max-members-per-community "$PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY" \
    --alpha "$PFM_HOMOGENEITY_ALPHA" \
    --outfile-prefix "ManualAdjustedConsensusHomogeneity" \
    --title "Manual-adjusted Infomap consensus homogeneity" \
    "${HOMOG_ROT_ARGS[@]}"
fi

echo "[pfm_update] complete"

#!/usr/bin/env bash
# PFM stage entrypoint: distance build, ridge fusion, and optional areal parcellation.
set -euo pipefail

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
MEDIR="${3:?missing MEDIR}"
_unused_start_session="${4:?missing StartSession}"
FuncDirName="${5:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${6:-${FUNC_FILE_PREFIX:-Rest}}"

SubjectDir="${StudyFolder}/${Subject}"
PFM_STRATEGY="${PFM_STRATEGY:-ridge_fusion}"
PFM_PYTHON="${PFM_PYTHON:-${PIPELINE_PYTHON:-python3}}"
PFM_RESOURCES_ROOT="${PFM_RESOURCES_ROOT:-${MEDIR}/res0urces}"
PFM_ROOT_DIR="${PFM_OUTDIR:-${SubjectDir}/func/${FuncDirName}/PFM}"
if [[ "$PFM_STRATEGY" == "ridge_fusion" ]]; then
  PFM_OUTDIR="${PFM_ROOT_DIR}/RidgeFusion"
else
  PFM_OUTDIR="${PFM_ROOT_DIR}/Infomap"
fi
PFM_PREP_DIR="${PFM_PREP_DIR:-}"
PFM_INPUT_CIFTI="${PFM_INPUT_CIFTI:-}"
PFM_INPUT_TAG="${PFM_INPUT_TAG:-${CONCAT_INPUT_TAG:-OCME+MEICA+MGTR}}"
PFM_CONCAT_OUT_SUBDIR="${PFM_CONCAT_OUT_SUBDIR:-${CONCAT_OUT_SUBDIR:-ConcatenatedCiftis}}"
PFM_FD_THRESHOLD="${PFM_FD_THRESHOLD:-${CONCAT_FD_THRESHOLD:-0.3}}"

PFM_DISTANCE_MATRIX="${PFM_DISTANCE_MATRIX:-${SubjectDir}/anat/T1w/fsaverage_LR32k/DistanceMatrix.npy}"
PFM_DISTANCE_BUILD_IF_MISSING="${PFM_DISTANCE_BUILD_IF_MISSING:-1}"
PFM_DISTANCE_VARIANT_CHUNK_ROWS="${PFM_DISTANCE_VARIANT_CHUNK_ROWS:-128}"
PFM_DISTANCE_CORTEX_MODE="${PFM_DISTANCE_CORTEX_MODE:-hybrid}"
PFM_DISTANCE_EUCLIDEAN_OVERRIDE_MM="${PFM_DISTANCE_EUCLIDEAN_OVERRIDE_MM:-5}"

PFM_RF_OUTFILE="${PFM_RF_OUTFILE:-RidgeFusion_VTX}"
PFM_RF_FC_WEIGHT="${PFM_RF_FC_WEIGHT:-1.0}"
PFM_RF_FC_DEMEAN="${PFM_RF_FC_DEMEAN:-0}"
PFM_RF_SPATIAL_WEIGHT="${PFM_RF_SPATIAL_WEIGHT:-0.1}"
PFM_RF_LAMBDA="${PFM_RF_LAMBDA:-10}"
PFM_RF_LOCAL_EXCLUSION_MM="${PFM_RF_LOCAL_EXCLUSION_MM:-10}"
PFM_RF_SUBCORT_REGRESS_ENABLE="${PFM_RF_SUBCORT_REGRESS_ENABLE:-1}"
PFM_RF_SUBCORT_REGRESS_DISTANCE_MM="${PFM_RF_SUBCORT_REGRESS_DISTANCE_MM:-20}"
PFM_RF_BRAIN_STRUCTURES_CSV="${PFM_RF_BRAIN_STRUCTURES_CSV:-CORTEX_LEFT,CEREBELLUM_LEFT,ACCUMBENS_LEFT,CAUDATE_LEFT,PUTAMEN_LEFT,THALAMUS_LEFT,HIPPOCAMPUS_LEFT,AMYGDALA_LEFT,CORTEX_RIGHT,CEREBELLUM_RIGHT,ACCUMBENS_RIGHT,CAUDATE_RIGHT,PUTAMEN_RIGHT,THALAMUS_RIGHT,HIPPOCAMPUS_RIGHT,AMYGDALA_RIGHT}"
PFM_RF_SMOOTHING_KERNEL="${PFM_RF_SMOOTHING_KERNEL:-1.7}"
PFM_PRIORS_MAT="${PFM_PRIORS_MAT:-${NETWORK_PRIORS_MAT:-}}"
PFM_SUBCORT_PRIORS_NII="${PFM_SUBCORT_PRIORS_NII:-}"
PFM_FUNC_XFMS_DIRNAME="${PFM_FUNC_XFMS_DIRNAME:-${FUNC_XFMS_DIRNAME:-${FuncDirName}}}"

PFM_AREAL_ENABLE="${PFM_AREAL_ENABLE:-0}"
PFM_AREAL_OUTFILE="${PFM_AREAL_OUTFILE:-}"
PFM_AREAL_MIN_SIZE="${PFM_AREAL_MIN_SIZE:-30}"
PFM_NEIGHBORS_MAT="${PFM_NEIGHBORS_MAT:-${PFM_RESOURCES_ROOT}/Cifti_surf_neighbors_LR_normalwall.mat}"

# Opt-in, experimental multi-network extension for Ridge-Fusion areal output.
# Preliminary behavior is encouraging, but broader validation is still needed.
PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE="${PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE:-0}"
PFM_RF_MULTINETWORK_STRIPE_PRESET="${PFM_RF_MULTINETWORK_STRIPE_PRESET:-balanced}"
PFM_RF_MULTINETWORK_OUTDIR_NAME="${PFM_RF_MULTINETWORK_OUTDIR_NAME:-ExperimentalMultiNetwork}"
PFM_RF_MULTINETWORK_OUTFILE="${PFM_RF_MULTINETWORK_OUTFILE:-${PFM_RF_OUTFILE}+ExperimentalMultiNetwork}"
PFM_RF_MULTINETWORK_REFERENCE_COUNT="${PFM_RF_MULTINETWORK_REFERENCE_COUNT:-1500}"
PFM_RF_MULTINETWORK_REFERENCE_SEED="${PFM_RF_MULTINETWORK_REFERENCE_SEED:-12345}"
PFM_RF_MULTINETWORK_CV_LOCAL_EXCLUSION_MM="${PFM_RF_MULTINETWORK_CV_LOCAL_EXCLUSION_MM:-30}"
PFM_RF_MULTINETWORK_SPLIT_REPEATS="${PFM_RF_MULTINETWORK_SPLIT_REPEATS:-8}"
PFM_RF_MULTINETWORK_SPLIT_BLOCK_LENGTH="${PFM_RF_MULTINETWORK_SPLIT_BLOCK_LENGTH:-25}"
PFM_RF_MULTINETWORK_SPLIT_SEED="${PFM_RF_MULTINETWORK_SPLIT_SEED:-12345}"
PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION="${PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION:-0.8}"
PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION="${PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION:-0.8}"
PFM_RF_MULTINETWORK_SELECTION_CORRECTION="${PFM_RF_MULTINETWORK_SELECTION_CORRECTION:-fdr}"
PFM_RF_MULTINETWORK_FDR_ALPHA="${PFM_RF_MULTINETWORK_FDR_ALPHA:-0.05}"
PFM_RF_MULTINETWORK_PARCEL_NULL_ITERATIONS="${PFM_RF_MULTINETWORK_PARCEL_NULL_ITERATIONS:-2000}"
PFM_RF_MULTINETWORK_PARCEL_NULL_SEED="${PFM_RF_MULTINETWORK_PARCEL_NULL_SEED:-12345}"
PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT="${PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT:-0.20}"
PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT="${PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT:-0.15}"
PFM_RF_MULTINETWORK_NULL_ALPHA="${PFM_RF_MULTINETWORK_NULL_ALPHA:-0.05}"
PFM_RF_MULTINETWORK_MAX_NETWORKS="${PFM_RF_MULTINETWORK_MAX_NETWORKS:-3}"
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN="${PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN:-0.10}"
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN="${PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN:-0.40}"
PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES="${PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES:-20}"
PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS="${PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS:-5}"

if [[ "$PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE" == "1" ]]; then
  case "$PFM_RF_MULTINETWORK_STRIPE_PRESET" in
    strict)
      PFM_RF_MULTINETWORK_SELECTION_CORRECTION="fwer"
      PFM_RF_MULTINETWORK_FDR_ALPHA="0.05"
      PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION="0.90"
      PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION="0.90"
      PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT="0.25"
      PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT="0.20"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN="0.20"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN="0.667"
      PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES="20"
      PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS="5"
      ;;
    balanced)
      PFM_RF_MULTINETWORK_SELECTION_CORRECTION="fdr"
      PFM_RF_MULTINETWORK_FDR_ALPHA="0.05"
      PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION="0.80"
      PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION="0.80"
      PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT="0.20"
      PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT="0.15"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN="0.10"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN="0.40"
      PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES="20"
      PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS="5"
      ;;
    loose)
      PFM_RF_MULTINETWORK_SELECTION_CORRECTION="fdr"
      PFM_RF_MULTINETWORK_FDR_ALPHA="0.10"
      PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION="0.70"
      PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION="0.70"
      PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT="0.15"
      PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT="0.10"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN="0.075"
      PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN="0.333"
      PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES="10"
      PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS="3"
      ;;
    custom)
      ;;
    *)
      echo "ERROR: PFM_RF_MULTINETWORK_STRIPE_PRESET must be strict, balanced, loose, or custom (got: $PFM_RF_MULTINETWORK_STRIPE_PRESET)"
      exit 2
      ;;
  esac
fi

if [[ -z "${PFM_INFOMAP_DISTANCE_MATRIX:-}" ]]; then
  PFM_INFOMAP_DISTANCE_MATRIX="$PFM_DISTANCE_MATRIX"
fi
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="${PFM_INFOMAP_GRAPH_DENSITIES_EXPR:-0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001}"
PFM_INFOMAP_NUM_REPS_EXPR="${PFM_INFOMAP_NUM_REPS_EXPR:-1,2,5,10,20,30,50,75,100}"
PFM_INFOMAP_MIN_DISTANCE="${PFM_INFOMAP_MIN_DISTANCE:-30}"
PFM_INFOMAP_BAD_VERTS_CSV="${PFM_INFOMAP_BAD_VERTS_CSV:-}"
PFM_INFOMAP_STRUCTURES_CSV="${PFM_INFOMAP_STRUCTURES_CSV:-}"
PFM_INFOMAP_NUM_CORES="${PFM_INFOMAP_NUM_CORES:-1}"
PFM_INFOMAP_BINARY="${PFM_INFOMAP_BINARY:-}"
PFM_INFOMAP_NETWORK_MAPPING_ENABLE="${PFM_INFOMAP_NETWORK_MAPPING_ENABLE:-0}"
PFM_INFOMAP_DRY_RUN="${PFM_INFOMAP_DRY_RUN:-0}"
PFM_INFOMAP_LABEL_OUTFILE="${PFM_INFOMAP_LABEL_OUTFILE:-InfomapNetworkLabels}"
PFM_INFOMAP_LABEL_FC_WEIGHT="${PFM_INFOMAP_LABEL_FC_WEIGHT:-1.0}"
PFM_INFOMAP_LABEL_SPATIAL_WEIGHT="${PFM_INFOMAP_LABEL_SPATIAL_WEIGHT:-1.0}"
PFM_INFOMAP_LABEL_CONFIDENCE_THRESHOLD="${PFM_INFOMAP_LABEL_CONFIDENCE_THRESHOLD:-0.15}"
PFM_INFOMAP_LABEL_MIN_FC_SIMILARITY="${PFM_INFOMAP_LABEL_MIN_FC_SIMILARITY:-0.33}"
PFM_INFOMAP_LABEL_MIN_COMMUNITY_SIZE="${PFM_INFOMAP_LABEL_MIN_COMMUNITY_SIZE:-10}"
PFM_INFOMAP_LABEL_UNASSIGNED_VALUE="${PFM_INFOMAP_LABEL_UNASSIGNED_VALUE:-21}"
PFM_INFOMAP_LABEL_STRICT_THRESHOLDING="${PFM_INFOMAP_LABEL_STRICT_THRESHOLDING:-0}"
PFM_INFOMAP_LABEL_DENSITY_INDEX="${PFM_INFOMAP_LABEL_DENSITY_INDEX:--1}"
PFM_INFOMAP_LABEL_WB_COMMAND="${PFM_INFOMAP_LABEL_WB_COMMAND:-wb_command}"
PFM_INFOMAP_LABEL_WRITE_FC="${PFM_INFOMAP_LABEL_WRITE_FC:-1}"
PFM_INFOMAP_LABEL_FC_CHUNK_ROWS="${PFM_INFOMAP_LABEL_FC_CHUNK_ROWS:-4096}"
PFM_INFOMAP_MANUAL_LABEL_APPLY_ENABLE="${PFM_INFOMAP_MANUAL_LABEL_APPLY_ENABLE:-0}"
PFM_INFOMAP_MANUAL_LABEL_TABLE="${PFM_INFOMAP_MANUAL_LABEL_TABLE:-}"
PFM_INFOMAP_MANUAL_LABEL_OUTFILE="${PFM_INFOMAP_MANUAL_LABEL_OUTFILE:-${PFM_INFOMAP_LABEL_OUTFILE}_ManualAdjusted}"
PFM_INFOMAP_UPDATE_ENABLE="${PFM_INFOMAP_UPDATE_ENABLE:-0}"
PFM_INFOMAP_UPDATE_TABLE_GLOB="${PFM_INFOMAP_UPDATE_TABLE_GLOB:-GraphDensity_*/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv}"
PFM_INFOMAP_UPDATE_OUTFILE="${PFM_INFOMAP_UPDATE_OUTFILE:-InfomapNetworkLabels_ManualAdjusted}"
PFM_HOMOGENEITY_TEST_ENABLE="${PFM_HOMOGENEITY_TEST_ENABLE:-0}"
PFM_HOMOGENEITY_ROTATIONS_CIFTI="${PFM_HOMOGENEITY_ROTATIONS_CIFTI:-${PFM_RESOURCES_ROOT}/Rotated_inds.dtseries.nii}"
PFM_HOMOGENEITY_N_ROTATIONS="${PFM_HOMOGENEITY_N_ROTATIONS:-1000}"
PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE="${PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE:-5}"
PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY="${PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY:-1000}"
PFM_HOMOGENEITY_ALPHA="${PFM_HOMOGENEITY_ALPHA:-0.05}"

if [[ "$PFM_STRATEGY" != "ridge_fusion" && "$PFM_STRATEGY" != "infomap" ]]; then
  echo "ERROR: PFM_STRATEGY must be ridge_fusion or infomap (got: $PFM_STRATEGY)"
  exit 2
fi

if [[ -z "$PFM_INPUT_CIFTI" ]]; then
  FD_TAG="${PFM_FD_THRESHOLD//./p}"
  PFM_INPUT_CIFTI="${SubjectDir}/func/${FuncDirName}/${PFM_CONCAT_OUT_SUBDIR}/${FuncFilePrefix}_${PFM_INPUT_TAG}_Concatenated+FDlt${FD_TAG}.dtseries.nii"
fi

L_MID="${SubjectDir}/anat/T1w/fsaverage_LR32k/${Subject}.L.midthickness.32k_fs_LR.surf.gii"
R_MID="${SubjectDir}/anat/T1w/fsaverage_LR32k/${Subject}.R.midthickness.32k_fs_LR.surf.gii"
L_SPHERE="${SubjectDir}/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.sphere.32k_fs_LR.surf.gii"
R_SPHERE="${SubjectDir}/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.sphere.32k_fs_LR.surf.gii"

echo "[pfm] strategy=${PFM_STRATEGY}"
echo "[pfm] input CIFTI=${PFM_INPUT_CIFTI}"
echo "[pfm] output dir=${PFM_OUTDIR}"

[[ -f "$PFM_INPUT_CIFTI" ]] || { echo "ERROR: missing input CIFTI: $PFM_INPUT_CIFTI"; exit 2; }
[[ -f "$L_MID" && -f "$R_MID" ]] || { echo "ERROR: missing midthickness surfaces"; exit 2; }
command -v wb_command >/dev/null 2>&1 || { echo "ERROR: wb_command not found"; exit 2; }
if [[ "$PFM_STRATEGY" == "ridge_fusion" && "$PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE" == "1" && "$PFM_AREAL_ENABLE" != "1" ]]; then
  echo "ERROR: PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE=1 requires PFM_AREAL_ENABLE=1"
  exit 2
fi

mkdir -p "$PFM_OUTDIR"
if [[ -z "$PFM_PREP_DIR" ]]; then
  PREP_DIR="${PFM_OUTDIR}/prep"
else
  PREP_DIR="$PFM_PREP_DIR"
fi
mkdir -p "$PREP_DIR"
echo "[pfm] prep dir=${PREP_DIR} (PFM-only intermediates)"

if [[ "$PFM_STRATEGY" == "ridge_fusion" ]]; then
  echo "[pfm] distance matrix=${PFM_DISTANCE_MATRIX}"
  if [[ ! -f "$PFM_DISTANCE_MATRIX" && "$PFM_DISTANCE_BUILD_IF_MISSING" == "1" ]]; then
    echo "[pfm] building distance matrix (default model) -> ${PFM_DISTANCE_MATRIX}"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_distance_matrix_build.py" \
      --ref-cifti "$PFM_INPUT_CIFTI" \
      --left-surf "$L_MID" \
      --right-surf "$R_MID" \
      --out-npy "$PFM_DISTANCE_MATRIX" \
      --chunk-rows "$PFM_DISTANCE_VARIANT_CHUNK_ROWS" \
      --cortex-distance-mode "$PFM_DISTANCE_CORTEX_MODE" \
      --euclidean-override-mm "$PFM_DISTANCE_EUCLIDEAN_OVERRIDE_MM"
  fi
  [[ -f "$PFM_DISTANCE_MATRIX" ]] || { echo "ERROR: distance matrix not found: $PFM_DISTANCE_MATRIX"; exit 2; }
fi

if [[ "${PFM_RF_SUBCORT_REGRESS_ENABLE}" == "1" ]]; then
  IN_BASENAME="$(basename "$PFM_INPUT_CIFTI" .dtseries.nii)"
  PFM_INPUT_CIFTI_REG="${PREP_DIR}/${IN_BASENAME}+SubcortRegression.dtseries.nii"
  echo "[pfm] running Python subcortical regression (distance=${PFM_RF_SUBCORT_REGRESS_DISTANCE_MM} mm) -> ${PFM_INPUT_CIFTI_REG}"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_subcort_regress.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --out-cifti "$PFM_INPUT_CIFTI_REG" \
    --left-surf "$L_MID" \
    --right-surf "$R_MID" \
    --distance-npy "$PFM_DISTANCE_MATRIX" \
    --distance-mm "$PFM_RF_SUBCORT_REGRESS_DISTANCE_MM"
  PFM_INPUT_CIFTI="$PFM_INPUT_CIFTI_REG"
fi

if awk "BEGIN{exit !(${PFM_RF_SMOOTHING_KERNEL} > 0)}"; then
  IN_BASENAME="$(basename "$PFM_INPUT_CIFTI" .dtseries.nii)"
  PFM_INPUT_CIFTI_SMOOTH="${PREP_DIR}/${IN_BASENAME}+SpatialSmoothing${PFM_RF_SMOOTHING_KERNEL}.dtseries.nii"
  echo "[pfm] smoothing input CIFTI with kernel=${PFM_RF_SMOOTHING_KERNEL} mm -> ${PFM_INPUT_CIFTI_SMOOTH}"
  wb_command -cifti-smoothing "$PFM_INPUT_CIFTI" "$PFM_RF_SMOOTHING_KERNEL" "$PFM_RF_SMOOTHING_KERNEL" COLUMN \
    "$PFM_INPUT_CIFTI_SMOOTH" -left-surface "$L_MID" -right-surface "$R_MID" -merged-volume
  PFM_INPUT_CIFTI="$PFM_INPUT_CIFTI_SMOOTH"
fi

echo "[pfm] prepared input CIFTI=${PFM_INPUT_CIFTI}"

PFM_NETWORK_DLABEL=""
PFM_AREAL_DISTANCE_MATRIX=""

if [[ "$PFM_STRATEGY" == "ridge_fusion" ]]; then

  [[ -f "$PFM_DISTANCE_MATRIX" ]] || { echo "ERROR: distance matrix not found: $PFM_DISTANCE_MATRIX"; exit 2; }
  [[ -f "$PFM_PRIORS_MAT" ]] || { echo "ERROR: missing PFM cortical network priors mat: $PFM_PRIORS_MAT"; exit 2; }

  PFM_SUBCORT_PRIORS_ACPC=""
  if [[ -n "$PFM_SUBCORT_PRIORS_NII" ]]; then
    [[ -f "$PFM_SUBCORT_PRIORS_NII" ]] || { echo "ERROR: missing PFM subcortical priors NIfTI: $PFM_SUBCORT_PRIORS_NII"; exit 2; }
    XFM_STANDARD2ACPC="${SubjectDir}/anat/MNINonLinear/xfms/standard2acpc_dc.nii.gz"
    ACPC_REF_FUNC="${SubjectDir}/func/xfms/${PFM_FUNC_XFMS_DIRNAME}/T1w_acpc_brain_func.nii.gz"
    [[ -f "$XFM_STANDARD2ACPC" ]] || { echo "ERROR: missing standard2acpc warp: $XFM_STANDARD2ACPC"; exit 2; }
    [[ -f "$ACPC_REF_FUNC" ]] || { echo "ERROR: missing ACPC functional reference volume: $ACPC_REF_FUNC"; exit 2; }
    PFM_SUBCORT_PRIORS_ACPC="${PREP_DIR}/SubcorticalPriors_acpc.nii.gz"
    echo "[pfm] warping subcortical priors MNI->ACPC -> ${PFM_SUBCORT_PRIORS_ACPC}"
    applywarp -i "$PFM_SUBCORT_PRIORS_NII" -o "$PFM_SUBCORT_PRIORS_ACPC" -r "$ACPC_REF_FUNC" -w "$XFM_STANDARD2ACPC"
  fi

  echo "[pfm] running Python ridge fusion"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_ridge_fusion.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --distance-npy "$PFM_DISTANCE_MATRIX" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --outdir "$PFM_OUTDIR" \
    --outfile "$PFM_RF_OUTFILE" \
    --fc-weight "$PFM_RF_FC_WEIGHT" \
    --fc-demean "$PFM_RF_FC_DEMEAN" \
    --spatial-weight "$PFM_RF_SPATIAL_WEIGHT" \
    --lambda "$PFM_RF_LAMBDA" \
    --local-exclusion-mm "$PFM_RF_LOCAL_EXCLUSION_MM" \
    --brain-structures-csv "$PFM_RF_BRAIN_STRUCTURES_CSV" \
    --subcort-priors-nii "$PFM_SUBCORT_PRIORS_ACPC" \
    --left-surf "$L_MID" \
    --right-surf "$R_MID"
  PFM_NETWORK_DLABEL="${PFM_OUTDIR}/${PFM_RF_OUTFILE}.dlabel.nii"
  PFM_AREAL_DISTANCE_MATRIX="$PFM_DISTANCE_MATRIX"
else
  if [[ "$PFM_INFOMAP_NETWORK_MAPPING_ENABLE" == "1" ]]; then
    [[ -f "$PFM_PRIORS_MAT" ]] || { echo "ERROR: missing PFM cortical network priors mat for infomap network mapping: $PFM_PRIORS_MAT"; exit 2; }
  fi
  echo "[pfm] infomap distance matrix=${PFM_INFOMAP_DISTANCE_MATRIX}"
  if [[ ! -f "$PFM_INFOMAP_DISTANCE_MATRIX" && "$PFM_DISTANCE_BUILD_IF_MISSING" == "1" ]]; then
    if [[ "$PFM_INFOMAP_DISTANCE_MATRIX" == *.npy ]]; then
      echo "[pfm] building infomap distance matrix (default model) -> ${PFM_INFOMAP_DISTANCE_MATRIX}"
      "$PFM_PYTHON" "$MEDIR/lib/pfm_distance_matrix_build.py" \
        --ref-cifti "$PFM_INPUT_CIFTI" \
        --left-surf "$L_MID" \
        --right-surf "$R_MID" \
        --out-npy "$PFM_INFOMAP_DISTANCE_MATRIX" \
        --chunk-rows "$PFM_DISTANCE_VARIANT_CHUNK_ROWS" \
        --cortex-distance-mode "$PFM_DISTANCE_CORTEX_MODE" \
        --euclidean-override-mm "$PFM_DISTANCE_EUCLIDEAN_OVERRIDE_MM"
    fi
  fi
  [[ -f "$PFM_INFOMAP_DISTANCE_MATRIX" ]] || { echo "ERROR: missing PFM_INFOMAP_DISTANCE_MATRIX: $PFM_INFOMAP_DISTANCE_MATRIX"; exit 2; }

  if [[ "$PFM_INFOMAP_NETWORK_MAPPING_ENABLE" == "1" ]]; then
    echo "[pfm] running Python infomap with community-to-network ID mapping"
  else
    echo "[pfm] running Python infomap (community mapping only)"
  fi
  INFOMAP_ARGS=(
    --in-cifti "$PFM_INPUT_CIFTI"
    --distance "$PFM_INFOMAP_DISTANCE_MATRIX"
    --outdir "$PFM_OUTDIR"
    --graph-densities "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR"
    --num-reps "$PFM_INFOMAP_NUM_REPS_EXPR"
    --min-distance "$PFM_INFOMAP_MIN_DISTANCE"
    --num-cores "$PFM_INFOMAP_NUM_CORES"
    --density-subdirs 1
  )
  if [[ -n "$PFM_INFOMAP_BINARY" ]]; then
    INFOMAP_ARGS+=( --infomap-binary "$PFM_INFOMAP_BINARY" )
  fi
  if [[ -n "$PFM_INFOMAP_STRUCTURES_CSV" ]]; then
    INFOMAP_ARGS+=( --structures-csv "$PFM_INFOMAP_STRUCTURES_CSV" )
  fi
  if [[ -n "$PFM_INFOMAP_BAD_VERTS_CSV" ]]; then
    INFOMAP_ARGS+=( --bad-verts-csv "$PFM_INFOMAP_BAD_VERTS_CSV" )
  fi
  if [[ "$PFM_INFOMAP_DRY_RUN" == "1" ]]; then
    INFOMAP_ARGS+=( --dry-run )
  fi
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    "$PFM_PYTHON" "$MEDIR/lib/pfm_infomap.py" "${INFOMAP_ARGS[@]}"

  if [[ "$PFM_INFOMAP_NETWORK_MAPPING_ENABLE" == "1" && "$PFM_INFOMAP_DRY_RUN" != "1" ]]; then
    echo "[pfm] labeling infomap communities with canonical network IDs"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_infomap_labeler.py" \
      --in-cifti "$PFM_INPUT_CIFTI" \
      --communities-cifti "${PFM_OUTDIR}/Bipartite_PhysicalCommunities.dtseries.nii" \
      --priors-mat "$PFM_PRIORS_MAT" \
      --outdir "$PFM_OUTDIR" \
      --outfile-prefix "$PFM_INFOMAP_LABEL_OUTFILE" \
      --density-index "$PFM_INFOMAP_LABEL_DENSITY_INDEX" \
      --density-values "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR" \
      --write-community-fc "$PFM_INFOMAP_LABEL_WRITE_FC" \
      --fc-chunk-rows "$PFM_INFOMAP_LABEL_FC_CHUNK_ROWS" \
      --fc-weight "$PFM_INFOMAP_LABEL_FC_WEIGHT" \
      --spatial-weight "$PFM_INFOMAP_LABEL_SPATIAL_WEIGHT" \
      --confidence-threshold "$PFM_INFOMAP_LABEL_CONFIDENCE_THRESHOLD" \
      --min-fc-similarity "$PFM_INFOMAP_LABEL_MIN_FC_SIMILARITY" \
      --min-community-size "$PFM_INFOMAP_LABEL_MIN_COMMUNITY_SIZE" \
      --unassigned-value "$PFM_INFOMAP_LABEL_UNASSIGNED_VALUE" \
      --strict-thresholding "$PFM_INFOMAP_LABEL_STRICT_THRESHOLDING" \
      --left-surf "$L_MID" \
      --right-surf "$R_MID" \
      --wb-command "$PFM_INFOMAP_LABEL_WB_COMMAND" \
      --density-output-mode subdirs
      PFM_NETWORK_DLABEL="${PFM_OUTDIR}/${PFM_INFOMAP_LABEL_OUTFILE}_ModeConsensus.dlabel.nii"

    if [[ "$PFM_INFOMAP_MANUAL_LABEL_APPLY_ENABLE" == "1" ]]; then
      if [[ -z "$PFM_INFOMAP_MANUAL_LABEL_TABLE" ]]; then
        PFM_INFOMAP_MANUAL_LABEL_TABLE=""
      fi
      echo "[pfm] applying manual infomap network label corrections"
      MANUAL_ARGS=()
      if [[ -n "$PFM_INFOMAP_MANUAL_LABEL_TABLE" ]]; then
        [[ -f "$PFM_INFOMAP_MANUAL_LABEL_TABLE" ]] || { echo "ERROR: missing Infomap manual correction table: $PFM_INFOMAP_MANUAL_LABEL_TABLE"; exit 2; }
        MANUAL_ARGS+=( --manual-corrections "$PFM_INFOMAP_MANUAL_LABEL_TABLE" )
      else
        MANUAL_ARGS+=( --manual-corrections-glob "$PFM_INFOMAP_UPDATE_TABLE_GLOB" )
      fi
      "$PFM_PYTHON" "$MEDIR/lib/pfm_infomap_manual_labels.py" \
        --communities-cifti "${PFM_OUTDIR}/Bipartite_PhysicalCommunities.dtseries.nii" \
        "${MANUAL_ARGS[@]}" \
        --priors-mat "$PFM_PRIORS_MAT" \
        --outdir "$PFM_OUTDIR" \
        --outfile-prefix "$PFM_INFOMAP_MANUAL_LABEL_OUTFILE" \
        --density-index "$PFM_INFOMAP_LABEL_DENSITY_INDEX" \
        --density-values "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR" \
        --unassigned-value "$PFM_INFOMAP_LABEL_UNASSIGNED_VALUE" \
        --left-surf "$L_MID" \
        --right-surf "$R_MID" \
        --wb-command "$PFM_INFOMAP_LABEL_WB_COMMAND" \
        --density-output-mode subdirs
      PFM_NETWORK_DLABEL="${PFM_OUTDIR}/${PFM_INFOMAP_MANUAL_LABEL_OUTFILE}_ModeConsensus.dlabel.nii"
    fi
  fi
  PFM_AREAL_DISTANCE_MATRIX="$PFM_INFOMAP_DISTANCE_MATRIX"
fi

if [[ "$PFM_HOMOGENEITY_TEST_ENABLE" == "1" ]]; then
  HOMOG_ROT_ARGS=()
  if [[ -f "$PFM_HOMOGENEITY_ROTATIONS_CIFTI" ]]; then
    HOMOG_ROT_ARGS=( --rotations-cifti "$PFM_HOMOGENEITY_ROTATIONS_CIFTI" --n-rotations "$PFM_HOMOGENEITY_N_ROTATIONS" )
  else
    echo "[pfm] homogeneity rotations not found; running observed-only: $PFM_HOMOGENEITY_ROTATIONS_CIFTI"
    HOMOG_ROT_ARGS=( --n-rotations 0 )
  fi
  if [[ "$PFM_STRATEGY" == "infomap" ]]; then
    echo "[pfm] running Infomap density homogeneity diagnostics"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_homogeneity_test.py" \
      --in-cifti "$PFM_INPUT_CIFTI" \
      --labels-cifti "${PFM_OUTDIR}/Bipartite_PhysicalCommunities.dtseries.nii" \
      --outdir "${PFM_OUTDIR}/HomogeneityTest" \
      --density-values "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR" \
      --min-community-size "$PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE" \
      --max-members-per-community "$PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY" \
      --alpha "$PFM_HOMOGENEITY_ALPHA" \
      --outfile-prefix "InfomapDensityHomogeneity" \
      --title "Infomap density homogeneity" \
      "${HOMOG_ROT_ARGS[@]}"
  elif [[ -n "$PFM_NETWORK_DLABEL" && -f "$PFM_NETWORK_DLABEL" ]]; then
    echo "[pfm] running Ridge-Fusion homogeneity diagnostics"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_homogeneity_test.py" \
      --in-cifti "$PFM_INPUT_CIFTI" \
      --labels-cifti "$PFM_NETWORK_DLABEL" \
      --outdir "${PFM_OUTDIR}/HomogeneityTest" \
      --min-community-size "$PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE" \
      --max-members-per-community "$PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY" \
      --alpha "$PFM_HOMOGENEITY_ALPHA" \
      --outfile-prefix "RidgeFusionHomogeneity" \
      --title "Ridge-Fusion homogeneity" \
      "${HOMOG_ROT_ARGS[@]}"
  fi
fi

if [[ "$PFM_AREAL_ENABLE" == "1" ]]; then
  [[ -n "$PFM_NETWORK_DLABEL" ]] || { echo "ERROR: PFM_AREAL_ENABLE=1 but no network-level dlabel was generated for strategy=${PFM_STRATEGY}"; exit 2; }
  [[ -f "$PFM_NETWORK_DLABEL" ]] || { echo "ERROR: missing network-level dlabel for areal parcellation: $PFM_NETWORK_DLABEL"; exit 2; }
  [[ -f "$PFM_AREAL_DISTANCE_MATRIX" ]] || { echo "ERROR: missing distance matrix for areal parcellation: $PFM_AREAL_DISTANCE_MATRIX"; exit 2; }
  [[ "$PFM_AREAL_DISTANCE_MATRIX" == *.npy ]] || { echo "ERROR: areal parcellation requires a .npy distance matrix, got: $PFM_AREAL_DISTANCE_MATRIX"; exit 2; }
  [[ -f "$PFM_PRIORS_MAT" ]] || { echo "ERROR: missing PFM cortical network priors mat for areal parcellation: $PFM_PRIORS_MAT"; exit 2; }

  if [[ -z "$PFM_AREAL_OUTFILE" ]]; then
    if [[ "$PFM_STRATEGY" == "ridge_fusion" ]]; then
      PFM_AREAL_OUTFILE="${PFM_RF_OUTFILE}+ArealParcellation"
    elif [[ "$PFM_INFOMAP_MANUAL_LABEL_APPLY_ENABLE" == "1" ]]; then
      PFM_AREAL_OUTFILE="${PFM_INFOMAP_MANUAL_LABEL_OUTFILE}+ArealParcellation"
    else
      PFM_AREAL_OUTFILE="${PFM_INFOMAP_LABEL_OUTFILE}+ArealParcellation"
    fi
  fi

  echo "[pfm] running Python areal parcellation on ${PFM_NETWORK_DLABEL}"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_areal_parcellation.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --wta-dlabel "$PFM_NETWORK_DLABEL" \
    --neighbors-mat "$PFM_NEIGHBORS_MAT" \
    --distance-npy "$PFM_AREAL_DISTANCE_MATRIX" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --outdir "$PFM_OUTDIR" \
    --outfile "$PFM_AREAL_OUTFILE" \
    --min-parcel-size "$PFM_AREAL_MIN_SIZE" \
    --left-surf "$L_MID" \
    --right-surf "$R_MID"
fi

if [[ "$PFM_STRATEGY" == "ridge_fusion" && "$PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE" == "1" ]]; then
  echo "[pfm] experimental multi-network stripe preset=${PFM_RF_MULTINETWORK_STRIPE_PRESET}"
  PFM_RF_MULTINETWORK_DIR="${PFM_OUTDIR}/${PFM_RF_MULTINETWORK_OUTDIR_NAME}"
  PFM_RF_MULTINETWORK_CV_PREFIX="${PFM_RF_MULTINETWORK_OUTFILE}CV"
  PFM_RF_MULTINETWORK_PARCELS="${PFM_OUTDIR}/${PFM_AREAL_OUTFILE}.dlabel.nii"
  PFM_RF_MULTINETWORK_PROB="${PFM_OUTDIR}/${PFM_RF_OUTFILE}_ProbMaps.dtseries.nii"
  PFM_RF_MULTINETWORK_FD="${SubjectDir}/func/${FuncDirName}/${PFM_CONCAT_OUT_SUBDIR}/FD.txt"
  PFM_RF_MULTINETWORK_SCAN_INFO="${SubjectDir}/func/${FuncDirName}/${PFM_CONCAT_OUT_SUBDIR}/ScanInfo.json"
  PFM_RF_MULTINETWORK_SUPPORT="${PFM_RF_MULTINETWORK_DIR}/${PFM_RF_MULTINETWORK_CV_PREFIX}_ParcelSummary.csv"

  [[ -f "$PFM_RF_MULTINETWORK_PARCELS" ]] || { echo "ERROR: missing areal dlabel for experimental multi-network analysis: $PFM_RF_MULTINETWORK_PARCELS"; exit 2; }
  [[ -f "$PFM_RF_MULTINETWORK_PROB" ]] || { echo "ERROR: missing Ridge-Fusion probability maps: $PFM_RF_MULTINETWORK_PROB"; exit 2; }
  [[ -f "$PFM_RF_MULTINETWORK_FD" ]] || { echo "ERROR: missing concatenation FD file: $PFM_RF_MULTINETWORK_FD"; exit 2; }
  [[ -f "$PFM_RF_MULTINETWORK_SCAN_INFO" ]] || { echo "ERROR: missing concatenation ScanInfo.json: $PFM_RF_MULTINETWORK_SCAN_INFO"; exit 2; }
  [[ -f "$L_SPHERE" && -f "$R_SPHERE" ]] || { echo "ERROR: missing 32k spherical surfaces for stripe rendering"; exit 2; }
  mkdir -p "$PFM_RF_MULTINETWORK_DIR"

  echo "[pfm] experimental Ridge-Fusion multi-network CV/null selection"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_parcel_run_split_cv.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --fd-txt "$PFM_RF_MULTINETWORK_FD" \
    --scan-info-json "$PFM_RF_MULTINETWORK_SCAN_INFO" \
    --fd-threshold "$PFM_FD_THRESHOLD" \
    --parcels-dlabel "$PFM_RF_MULTINETWORK_PARCELS" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --distance-npy "$PFM_AREAL_DISTANCE_MATRIX" \
    --reference-count "$PFM_RF_MULTINETWORK_REFERENCE_COUNT" \
    --reference-seed "$PFM_RF_MULTINETWORK_REFERENCE_SEED" \
    --local-exclusion-mm "$PFM_RF_MULTINETWORK_CV_LOCAL_EXCLUSION_MM" \
    --split-repeats "$PFM_RF_MULTINETWORK_SPLIT_REPEATS" \
    --split-block-length "$PFM_RF_MULTINETWORK_SPLIT_BLOCK_LENGTH" \
    --split-seed "$PFM_RF_MULTINETWORK_SPLIT_SEED" \
    --min-consistency-fraction "$PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION" \
    --min-positive-gain-fraction "$PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION" \
    --selection-correction "$PFM_RF_MULTINETWORK_SELECTION_CORRECTION" \
    --fdr-alpha "$PFM_RF_MULTINETWORK_FDR_ALPHA" \
    --parcel-null-iterations "$PFM_RF_MULTINETWORK_PARCEL_NULL_ITERATIONS" \
    --parcel-null-seed "$PFM_RF_MULTINETWORK_PARCEL_NULL_SEED" \
    --min-second-network-weight "$PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT" \
    --min-third-network-weight "$PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT" \
    --null-alpha "$PFM_RF_MULTINETWORK_NULL_ALPHA" \
    --outdir "$PFM_RF_MULTINETWORK_DIR" \
    --outfile-prefix "$PFM_RF_MULTINETWORK_CV_PREFIX"

  echo "[pfm] rendering null-selected whole-parcel network stripes"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_parcel_uncertainty_stripes.py" \
    --parcels-dlabel "$PFM_RF_MULTINETWORK_PARCELS" \
    --prob-maps-cifti "$PFM_RF_MULTINETWORK_PROB" \
    --parcel-support-csv "$PFM_RF_MULTINETWORK_SUPPORT" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --left-sphere "$L_SPHERE" \
    --right-sphere "$R_SPHERE" \
    --max-networks-per-parcel "$PFM_RF_MULTINETWORK_MAX_NETWORKS" \
    --min-parcel-vertices "$PFM_AREAL_MIN_SIZE" \
    --min-subcortical-parcel-voxels "$PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS" \
    --outdir "$PFM_RF_MULTINETWORK_DIR" \
    --outfile-prefix "${PFM_RF_MULTINETWORK_OUTFILE}_Parcel"

  echo "[pfm] localizing null-selected support to contiguous vertex/voxel clusters"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_vertex_uncertainty_stripes.py" \
    --mode parcel_gated \
    --wta-dlabel "$PFM_NETWORK_DLABEL" \
    --prob-maps-cifti "$PFM_RF_MULTINETWORK_PROB" \
    --parcels-dlabel "$PFM_RF_MULTINETWORK_PARCELS" \
    --parcel-support-csv "$PFM_RF_MULTINETWORK_SUPPORT" \
    --neighbors-mat "$PFM_NEIGHBORS_MAT" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --left-sphere "$L_SPHERE" \
    --right-sphere "$R_SPHERE" \
    --secondary-prob-min "$PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN" \
    --secondary-to-top-ratio-min "$PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN" \
    --max-networks-per-zone "$PFM_RF_MULTINETWORK_MAX_NETWORKS" \
    --min-cortical-cluster-vertices "$PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES" \
    --min-subcortical-cluster-voxels "$PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS" \
    --outdir "$PFM_RF_MULTINETWORK_DIR" \
    --outfile-prefix "${PFM_RF_MULTINETWORK_OUTFILE}_Localized"
fi

echo "[pfm] complete"

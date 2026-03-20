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
PFM_PYTHON="${PFM_PYTHON:-python3}"
PFM_RESOURCES_ROOT="${PFM_RESOURCES_ROOT:-${MEDIR}/res0urces}"
PFM_INFOMAP_WRAPPER="${PFM_INFOMAP_WRAPPER:-${PFM_RESOURCES_ROOT}/PFM-InfoMap-Tmp/pfm_wrapper.m}"
PFM_OUTDIR="${PFM_OUTDIR:-${SubjectDir}/func/${FuncDirName}/PFM}"
PFM_INPUT_CIFTI="${PFM_INPUT_CIFTI:-}"
PFM_INPUT_TAG="${PFM_INPUT_TAG:-${CONCAT_INPUT_TAG:-OCME+MEICA+MGTR}}"
PFM_CONCAT_OUT_SUBDIR="${PFM_CONCAT_OUT_SUBDIR:-${CONCAT_OUT_SUBDIR:-ConcatenatedCiftis}}"
PFM_FD_THRESHOLD="${PFM_FD_THRESHOLD:-${CONCAT_FD_THRESHOLD:-0.3}}"

PFM_DISTANCE_MATRIX="${PFM_DISTANCE_MATRIX:-${SubjectDir}/anat/T1w/fsaverage_LR32k/DistanceMatrix.npy}"
PFM_DISTANCE_BUILD_IF_MISSING="${PFM_DISTANCE_BUILD_IF_MISSING:-1}"
PFM_DISTANCE_VARIANT_CHUNK_ROWS="${PFM_DISTANCE_VARIANT_CHUNK_ROWS:-128}"

PFM_RF_OUTFILE="${PFM_RF_OUTFILE:-RidgeFusion_VTX}"
PFM_RF_FC_WEIGHT="${PFM_RF_FC_WEIGHT:-1.0}"
PFM_RF_SPATIAL_WEIGHT="${PFM_RF_SPATIAL_WEIGHT:-0.1}"
PFM_RF_LAMBDA="${PFM_RF_LAMBDA:-10}"
PFM_RF_LOCAL_EXCLUSION_MM="${PFM_RF_LOCAL_EXCLUSION_MM:-10}"
PFM_RF_SUBCORT_REGRESS_ENABLE="${PFM_RF_SUBCORT_REGRESS_ENABLE:-1}"
PFM_RF_SUBCORT_REGRESS_DISTANCE_MM="${PFM_RF_SUBCORT_REGRESS_DISTANCE_MM:-20}"
PFM_RF_BRAIN_STRUCTURES="${PFM_RF_BRAIN_STRUCTURES:-}"
PFM_RF_BRAIN_STRUCTURES_CSV="${PFM_RF_BRAIN_STRUCTURES_CSV:-CORTEX_LEFT,CEREBELLUM_LEFT,ACCUMBENS_LEFT,CAUDATE_LEFT,PUTAMEN_LEFT,THALAMUS_LEFT,HIPPOCAMPUS_LEFT,AMYGDALA_LEFT,CORTEX_RIGHT,CEREBELLUM_RIGHT,ACCUMBENS_RIGHT,CAUDATE_RIGHT,PUTAMEN_RIGHT,THALAMUS_RIGHT,HIPPOCAMPUS_RIGHT,AMYGDALA_RIGHT}"
PFM_RF_SMOOTHING_KERNEL="${PFM_RF_SMOOTHING_KERNEL:-0}"
PFM_PRIORS_MAT="${PFM_PRIORS_MAT:-${NETWORK_PRIORS_MAT:-}}"

PFM_AREAL_ENABLE="${PFM_AREAL_ENABLE:-0}"
PFM_AREAL_OUTFILE="${PFM_AREAL_OUTFILE:-RidgeFusion_VTX+ArealParcellation}"
PFM_AREAL_MIN_SIZE="${PFM_AREAL_MIN_SIZE:-10}"
PFM_NEIGHBORS_MAT="${PFM_NEIGHBORS_MAT:-${PFM_RESOURCES_ROOT}/Cifti_surf_neighbors_LR_normalwall.mat}"

PFM_INFOMAP_MATLAB="${PFM_INFOMAP_MATLAB:-matlab}"
PFM_INFOMAP_DISTANCE_MATRIX="${PFM_INFOMAP_DISTANCE_MATRIX:-${SubjectDir}/anat/T1w/fsaverage_LR32k/DistanceMatrix.mat}"
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="${PFM_INFOMAP_GRAPH_DENSITIES_EXPR:-0.001:0.001:0.009}"
PFM_INFOMAP_NUM_REPS_EXPR="${PFM_INFOMAP_NUM_REPS_EXPR:-20}"
PFM_INFOMAP_MIN_DISTANCE="${PFM_INFOMAP_MIN_DISTANCE:-20}"
PFM_INFOMAP_BAD_VERTS_EXPR="${PFM_INFOMAP_BAD_VERTS_EXPR:-[]}"
PFM_INFOMAP_STRUCTURES_EXPR="${PFM_INFOMAP_STRUCTURES_EXPR:-[]}"
PFM_INFOMAP_BAD_VERTS_CSV="${PFM_INFOMAP_BAD_VERTS_CSV:-}"
PFM_INFOMAP_STRUCTURES_CSV="${PFM_INFOMAP_STRUCTURES_CSV:-}"
PFM_INFOMAP_NUM_CORES="${PFM_INFOMAP_NUM_CORES:-1}"
PFM_INFOMAP_BINARY="${PFM_INFOMAP_BINARY:-}"
PFM_INFOMAP_NETWORK_MAPPING_ENABLE="${PFM_INFOMAP_NETWORK_MAPPING_ENABLE:-0}"
PFM_INFOMAP_USE_MATLAB="${PFM_INFOMAP_USE_MATLAB:-0}"
PFM_INFOMAP_DRY_RUN="${PFM_INFOMAP_DRY_RUN:-0}"

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

echo "[pfm] strategy=${PFM_STRATEGY}"
echo "[pfm] input CIFTI=${PFM_INPUT_CIFTI}"
echo "[pfm] output dir=${PFM_OUTDIR}"

[[ -f "$PFM_INPUT_CIFTI" ]] || { echo "ERROR: missing input CIFTI: $PFM_INPUT_CIFTI"; exit 2; }
[[ -f "$L_MID" && -f "$R_MID" ]] || { echo "ERROR: missing midthickness surfaces"; exit 2; }
command -v wb_command >/dev/null 2>&1 || { echo "ERROR: wb_command not found"; exit 2; }

mkdir -p "$PFM_OUTDIR"
PREP_DIR="${PFM_OUTDIR}/prep"
mkdir -p "$PREP_DIR"

if [[ "$PFM_STRATEGY" == "ridge_fusion" ]]; then
  if [[ "${PFM_RF_SUBCORT_REGRESS_ENABLE}" == "1" ]]; then
    IN_BASENAME="$(basename "$PFM_INPUT_CIFTI" .dtseries.nii)"
    PFM_INPUT_CIFTI_REG="${PREP_DIR}/${IN_BASENAME}+SubcortRegression.dtseries.nii"
    echo "[pfm] running Python subcortical regression (distance=${PFM_RF_SUBCORT_REGRESS_DISTANCE_MM} mm) -> ${PFM_INPUT_CIFTI_REG}"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_subcort_regress.py" \
      --in-cifti "$PFM_INPUT_CIFTI" \
      --out-cifti "$PFM_INPUT_CIFTI_REG" \
      --left-surf "$L_MID" \
      --right-surf "$R_MID" \
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

  echo "[pfm] distance matrix=${PFM_DISTANCE_MATRIX}"
  if [[ ! -f "$PFM_DISTANCE_MATRIX" && "$PFM_DISTANCE_BUILD_IF_MISSING" == "1" ]]; then
    echo "[pfm] building distance matrix (default model) -> ${PFM_DISTANCE_MATRIX}"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_distance_matrix_build.py" \
      --ref-cifti "$PFM_INPUT_CIFTI" \
      --left-surf "$L_MID" \
      --right-surf "$R_MID" \
      --out-npy "$PFM_DISTANCE_MATRIX" \
      --chunk-rows "$PFM_DISTANCE_VARIANT_CHUNK_ROWS"
  fi
  [[ -f "$PFM_DISTANCE_MATRIX" ]] || { echo "ERROR: distance matrix not found: $PFM_DISTANCE_MATRIX"; exit 2; }
  [[ -f "$PFM_PRIORS_MAT" ]] || { echo "ERROR: missing PFM cortical network priors mat: $PFM_PRIORS_MAT"; exit 2; }

  echo "[pfm] running Python ridge fusion"
  "$PFM_PYTHON" "$MEDIR/lib/pfm_ridge_fusion.py" \
    --in-cifti "$PFM_INPUT_CIFTI" \
    --distance-npy "$PFM_DISTANCE_MATRIX" \
    --priors-mat "$PFM_PRIORS_MAT" \
    --outdir "$PFM_OUTDIR" \
    --outfile "$PFM_RF_OUTFILE" \
    --fc-weight "$PFM_RF_FC_WEIGHT" \
    --spatial-weight "$PFM_RF_SPATIAL_WEIGHT" \
    --lambda "$PFM_RF_LAMBDA" \
    --local-exclusion-mm "$PFM_RF_LOCAL_EXCLUSION_MM" \
    --brain-structures-csv "$PFM_RF_BRAIN_STRUCTURES_CSV" \
    --left-surf "$L_MID" \
    --right-surf "$R_MID"

  if [[ "$PFM_AREAL_ENABLE" == "1" ]]; then
    echo "[pfm] running Python areal parcellation"
    "$PFM_PYTHON" "$MEDIR/lib/pfm_areal_parcellation.py" \
      --in-cifti "$PFM_INPUT_CIFTI" \
      --wta-dlabel "${PFM_OUTDIR}/${PFM_RF_OUTFILE}.dlabel.nii" \
      --neighbors-mat "$PFM_NEIGHBORS_MAT" \
      --outdir "$PFM_OUTDIR" \
      --outfile "$PFM_AREAL_OUTFILE" \
      --min-size "$PFM_AREAL_MIN_SIZE" \
      --left-surf "$L_MID" \
      --right-surf "$R_MID"
  fi
else
  [[ "$PFM_INFOMAP_NETWORK_MAPPING_ENABLE" == "0" ]] || {
    echo "ERROR: PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1 is not supported yet."
    echo "       Current infomap path stops at community mapping outputs."
    exit 2
  }
  [[ -f "$PFM_INFOMAP_DISTANCE_MATRIX" ]] || { echo "ERROR: missing PFM_INFOMAP_DISTANCE_MATRIX: $PFM_INFOMAP_DISTANCE_MATRIX"; exit 2; }
  if [[ "${PFM_INFOMAP_USE_MATLAB}" == "1" ]]; then
    [[ -f "$PFM_INFOMAP_WRAPPER" ]] || { echo "ERROR: missing PFM_INFOMAP_WRAPPER: $PFM_INFOMAP_WRAPPER"; exit 2; }
    command -v "$PFM_INFOMAP_MATLAB" >/dev/null 2>&1 || { echo "ERROR: ${PFM_INFOMAP_MATLAB} not found"; exit 2; }

    WRAPPER_DIR="$(dirname "$PFM_INFOMAP_WRAPPER")"
    INFOMAP_BIN_EXPR="[]"
    if [[ -n "$PFM_INFOMAP_BINARY" ]]; then
      INFOMAP_BIN_EXPR="'${PFM_INFOMAP_BINARY}'"
    fi

    MAT_CMD="addpath(genpath('${PFM_RESOURCES_ROOT}/Utilities/msc/fieldtrip')); addpath('${WRAPPER_DIR}'); C=ft_read_cifti_mod('${PFM_INPUT_CIFTI}'); D=smartload('${PFM_INFOMAP_DISTANCE_MATRIX}'); pfm_wrapper(C,D,'${PFM_OUTDIR}',${PFM_INFOMAP_GRAPH_DENSITIES_EXPR},${PFM_INFOMAP_NUM_REPS_EXPR},${PFM_INFOMAP_MIN_DISTANCE},${PFM_INFOMAP_BAD_VERTS_EXPR},${PFM_INFOMAP_STRUCTURES_EXPR},${PFM_INFOMAP_NUM_CORES},${INFOMAP_BIN_EXPR});"
    echo "[pfm] running MATLAB infomap wrapper fallback (community mapping only)"
    "$PFM_INFOMAP_MATLAB" -batch "$MAT_CMD"
  else
    echo "[pfm] running Python infomap (community mapping only; no network identity mapping)"
    INFOMAP_ARGS=(
      --in-cifti "$PFM_INPUT_CIFTI"
      --distance "$PFM_INFOMAP_DISTANCE_MATRIX"
      --outdir "$PFM_OUTDIR"
      --graph-densities "$PFM_INFOMAP_GRAPH_DENSITIES_EXPR"
      --num-reps "$PFM_INFOMAP_NUM_REPS_EXPR"
      --min-distance "$PFM_INFOMAP_MIN_DISTANCE"
      --num-cores "$PFM_INFOMAP_NUM_CORES"
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
  fi
fi

echo "[pfm] complete"

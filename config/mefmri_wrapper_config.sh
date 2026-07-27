#!/usr/bin/env bash
# Config for ME-fMRI pipeline
#
# Wrapper usage:
#   mefmri_pipeline.sh <SubjectDir> [ConfigFile]
#
# This file is the main pipeline configuration. Most users should review the
# path/environment settings near the top, then change only task naming,
# processing mode, denoising method, and output-stage controls as needed.
#
# New user TODO before running:
#   - Confirm EnvironmentScript points to your HCP Pipelines setup script.
#   - Set CHARM_BIN to your SimNIBS CHARM executable.
#   - Confirm repository resources are present: AtlasTemplate, Priors.mat,
#     SubcorticalPriors.nii.gz, and template/resource files under res0urces/.
#   - Confirm conda environments or executable paths for tedana/ME-ICA,
#     ICA-AROMA, and experimental MEDIC/warpkit if those branches are used.
#   - Replace explicit tool overrides such as AROMA_BIN_OVERRIDE when using
#     ICA-AROMA.
#
# Organization:
#   1) Core paths and environment
#   2) Commonly used global knobs
#   3) Module-specific knobs
#   4) Advanced / rarely changed knobs

# =============================================================================
# 1) Core Paths and Environment
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EnvironmentScript="$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh"
CHARM_BIN="/home/charleslynch/SimNIBS-4.5/bin/charm"

# Repository resources used by downstream modules. Leave these repository-
# relative defaults unless your resource bundle lives somewhere else.
PFM_NETWORK_PRIORS_CORTICAL_MAT="$MEDIR/res0urces/Priors.mat"
PFM_NETWORK_PRIORS_SUBCORTICAL_NII="$MEDIR/res0urces/SubcorticalPriors.nii.gz"
PIPELINE_PYTHON="python3"

# =============================================================================
# 2) Commonly Used Global Knobs
# =============================================================================
# Resume / routing. Leave START_FROM_MODULE at validate for a fresh run.
START_SESSION=1
START_FROM_MODULE="validate"   # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|denoise|mgtr|vol2surf|concat|nsi|pfm|pfm_update
STOP_AFTER_MODULE=""           # "" to run full chain

# Functional naming and reference space
FUNC_DIRNAME="rest"              # task folder under func/unprocessed/; use "All" to loop over discovered tasks
FUNC_FILE_PREFIX="Rest"
FUNC_XFMS_DIRNAME=""             # empty => follow FUNC_DIRNAME; set explicitly only if you need a shared xfm namespace
DOF=6
AtlasTemplate="$MEDIR/res0urces/MNI152_T1_2mm.nii.gz"
AtlasSpace="T1w"                 # T1w|MNINonlinear
APPLY_N4_BIAS=0                  # 0|1 ; usually leave off if prescan normalize was enabled

# Global pipeline behavior. PROCESSING_MODE=auto chooses single_echo when fewer
# than three echoes are present, otherwise multi_echo.
RUN_CONFIG_SNAPSHOT=1            # 1 writes effective run metadata snapshot into func/<FUNC_DIRNAME>/qa
VALIDATE_ECHO_DIM4_POLICY="error" # error|warn
FUNC_NOFIELDMAP_MODE=0           # 0|1: functional-only fallback mode with zero-unwarp placeholders
PROCESSING_MODE="auto"          # auto|multi_echo|single_echo
MULTI_ECHO_DENOISE_METHOD="meica" # meica|aroma|acompcor
SINGLE_ECHO_DENOISE_METHOD="acompcor" # aroma|acompcor
SINGLE_ECHO_ECHO_INDEX=1        # fallback/source echo used for single-echo denoising
AROMA_NSI_THRESHOLD=0.05        # fixed NSI kill threshold for ICA-AROMA component screening
CONCAT_ENABLE=1
NSI_ENABLE=1                     # optional downstream metric; set 0 to skip when only preprocessing/QC is needed
PFM_ENABLE=1                     # optional downstream mapping; set 0 to skip time-consuming PFM outputs

# Threading defaults used by modules that support parallelism.
THREADS_DEFAULT=8
THREADS_ANAT_HCP="$THREADS_DEFAULT"
THREADS_FIELDMAPS="$THREADS_DEFAULT"
THREADS_COREG="$THREADS_DEFAULT"
THREADS_HEADMOTION="$THREADS_DEFAULT"
THREADS_MEICA=2

# =============================================================================
# 3) Module-Specific Knobs
# =============================================================================

# -----------------------------------------------------------------------------
# 3a) Anatomy / CHARM / HCP
# -----------------------------------------------------------------------------
HCP_ANAT_CLEAN_START=1           # 1 removes prior anat outputs before rerun
HCP_REGNAME="MSMSulc"            # MSMSulc|FS ; set FS if MSM isn't available
CHARM_REUSE_EXISTING_M2M=0       # 1 reuses existing anat/m2m_<Subject> outputs if present

# CHARM anatomical mask mode used to regenerate T1/T2 brain images.
# - charm: CHARM labeling-derived mask (<100), then dilated [default]
# - hcp: preserve pre-CHARM HCP-style whole-brain mask
CHARM_BRAIN_MASK_MODE="charm"    # charm|hcp
CHARM_BRAIN_MASK_DILATE_ITERS=1  # integer >= 0

# Cortical ribbon generation in anat module
CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS=1  # 1 removes HC/Amy/CSF, 0 keeps FS ribbon-only labels
CHARM_WRITE_CORTICAL_RIBBON=1           # 1 writes anat/T1w/CorticalRibbon.nii.gz

# -----------------------------------------------------------------------------
# 3b) Functional Distortion Correction / Coreg / Headmotion
# -----------------------------------------------------------------------------
DISTORTION_CORRECTION_MODE="topup" # topup|direct_b0|medic|none ; MEDIC/warpkit support is experimental
FM_PE_MODE="infer"               # infer|config
FM_AP_PE_DIR=""                  # e.g. j-
FM_PA_PE_DIR=""                  # e.g. j
TOPUP_CONFIG="b02b0.cnf"         # topup config name/path; for odd dimensions use "$MEDIR/templates/FSL/b02b0_nosubsamp.cnf"
EPIREG_PEDIR=""                  # empty => infer from PE.txt
SCAN_SPECIFIC_FM=1
FM_DIRECT_B0_DELTA_TE_MS=2.65    # fsl_prepare_fieldmap delta TE for direct B0/gradient-echo maps
FM_DIRECT_B0_PHASE_PATTERN="FieldMap_S*_R*.nii.gz,*_phasediff.nii.gz"
FM_DIRECT_B0_MAG_PATTERN="Mag_S*_R*.nii.gz,*_magnitude1.nii.gz"
HEADMOTION_KEEP_MCF=0
SBREF_FALLBACK_SKIP_TRS=10       # used only when no real SBref exists; also written to rmVols.txt
SBREF_REORIENT_TO_STD=1          # 1 runs fslreorient2std on real/fallback SBrefs before coreg
SBREF_ECHO_COMBINATION="mean"    # mean|sos|t2s|paid|first|last; used by topup, medic, and no-fieldmap coreg SBRefs
SBREF_MAX_TE_MS=60               # echo SBRefs with TE < this value are included; set to none to include all echoes
SBREF_T2SMAP_MASK_THR=1          # temporary nonzero mask threshold for t2s/paid SBRef combination
SBREF_T2SMAP_THREADS=1           # tedana t2smap threads for t2s/paid SBRef combination
SBREF_KEEP_INTERMEDIATES=0       # 1 keeps SBRef combination scratch folders for debugging
FUNC_REORIENT_TO_STD=1           # 1 runs fslreorient2std on copied echo timeseries before mcflirt

# Experimental MEDIC/warpkit runtime. Keep this separate from TEDANA_ENV so
# tedana's tested Python stack remains stable while warpkit can use Python >=3.11.
WARPKIT_ENV="warpkit_env"         # conda environment expected to provide wk-medic and related wk-* CLIs
WARPKIT_ACTIVATE_MODE="conda_activate" # conda_activate|conda_run|direct
WARPKIT_MEDIC_BIN="wk-medic"
WARPKIT_APPLY_WARP_BIN="wk-apply-warp"
WARPKIT_COMPUTE_JACOBIAN_BIN="wk-compute-jacobian"

# -----------------------------------------------------------------------------
# 3c) MGTR / Vol2Surf
# -----------------------------------------------------------------------------
# MGTR cortical ribbon source
# - xfms: func/xfms/<FUNC_XFMS_DIRNAME or FUNC_DIRNAME>/CorticalRibbon_*_func_mask.nii.gz
# - legacy_rois: func/rois/CorticalRibbon.nii.gz built from FreeSurfer ribbon mgz
MGTR_RIBBON_SOURCE="xfms"        # xfms|legacy_rois

MGTR_ENABLE="auto"               # auto|0|1; auto skips MGTR after aCompCor because ribbon GS is already modeled
MGTR_INPUT_TAG=""                # empty => derive from effective denoising branch
MGTR_OUTPUT_TAG=""               # empty => <MGTR_INPUT_TAG>+MGTR
VOL2SURF_INPUTS=""               # empty => derive from effective denoising branch
VOL2SURF_USE_CORTICAL_RIBBON_MASK=1
VOL2SURF_CIFTI_STAMP=""          # optional non-canonical suffix for experimental comparisons

# -----------------------------------------------------------------------------
# 3d) Tedana / MEICA / Reclassify
# -----------------------------------------------------------------------------
TEDANA_ENV="mefmri_env"          # conda environment used for tedana/ME-ICA
TEDANA_ACTIVATE_MODE="conda_activate" # conda_activate|conda_run|direct
TEDANA_COMPAT_MODE="modern"
MEPCA="350"
MaxIterations=500
MaxRestarts=5
TEDANA_VERBOSE=0
MEICA_SKIP_TEDANA_IF_EXISTS=0

MEICA_PARALLEL_JOBS=4
MEICA_RECLASSIFY_ENABLE=1
MEICA_TEDANA_SUBDIR="Tedana"
MEICA_RECLASS_SUBDIR="Reclassify"
MEICA_CLASSIFIER_MODE="nsi"      # nsi|legacy_template_rho|none
MEICA_RECLASS_NO_REPORTS=0       # 0 writes reports, 1 suppresses report generation
MEICA_QC_CIFTI_ENABLE=1
MEICA_QC_CIFTI_TAGS="betas_OC,t2sv,s0v"
MEICA_ORIG_ALIAS_ENABLE=1        # 1 creates legacy/orig filename aliases

# -----------------------------------------------------------------------------
# 3e) ICA-AROMA
# -----------------------------------------------------------------------------
AROMA_ENV="aroma"                # conda environment used for ICA-AROMA
AROMA_BIN_OVERRIDE="/home/charleslynch/aroma/ICA_AROMA.py"
AROMA_PYTHON="python"
AROMA_ACTIVATE_MODE="conda_activate" # conda_activate|conda_run|direct
AROMA_PARALLEL_JOBS=4
AROMA_OUT_SUBDIR="Aroma"
AROMA_OVERWRITE=1                # 0|1
AROMA_PLOT_REPORTS=0             # 0|1
AROMA_FEATURES_OUTPUT=1          # 0|1
AROMA_CLEAN_TYPE="nonaggr"       # nonaggr|aggr

# -----------------------------------------------------------------------------
# 3f) Concat
# -----------------------------------------------------------------------------
CONCAT_PYTHON="$PIPELINE_PYTHON"
CONCAT_INPUT_TAG=""              # empty => derive from effective denoising branch
CONCAT_OUT_SUBDIR="ConcatenatedCiftis"
CONCAT_CENSOR_BY_FD=1
CONCAT_FD_THRESHOLD=0.3

# Optional concat-stage QC movie. This runs immediately after concat when
# CRAWLING_SEED_FC_ENABLE=1; it is not a START/STOP routing module.
# It is resource-heavy: it creates a dense correlation CIFTI that can be tens of
# GB, then renders frames with Workbench/wb_surfer2. Expect tens of minutes to
# 1+ hour per subject; set CRAWLING_SEED_FC_ENABLE=0 to skip it.
CRAWLING_SEED_FC_ENABLE=0
CRAWLING_SEED_FC_INPUT_CIFTI=""          # empty => final concatenated CIFTI
CRAWLING_SEED_FC_INPUT_TAG=""            # empty => CONCAT_INPUT_TAG
CRAWLING_SEED_FC_OUTDIR=""               # empty => func/<FUNC_DIRNAME>/qa/CrawlingSeedFC
CRAWLING_SEED_FC_WB_SURFER2=""           # direct path to wb_surfer2; overrides conda env
CRAWLING_SEED_FC_WB_SURFER2_CONDA_ENV="wbsurfer_env"
CRAWLING_SEED_FC_SCENE_TEMPLATE="$MEDIR/res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene"
CRAWLING_SEED_FC_VERTICES="$MEDIR/res0urces/CrawlingSeedFC/VerticesToSample.txt"
CRAWLING_SEED_FC_TEMPLATE_SUBJECT="sub-TEMPLATE"
CRAWLING_SEED_FC_SURFACE_RESOURCE_DIR="" # empty => anat/T1w or anat/MNINonLinear fsaverage_LR32k from AtlasSpace
CRAWLING_SEED_FC_SURFACE_SUBJECT_PREFIX="" # empty => infer from subject surface filenames
CRAWLING_SEED_FC_FLAT_SURFACE_RESOURCE_DIR="" # empty => anat/MNINonLinear/fsaverage_LR32k for T1w flatmaps
CRAWLING_SEED_FC_FLAT_SURFACE_SUBJECT_PREFIX="" # empty => subject ID
CRAWLING_SEED_FC_WIDTH=1280
CRAWLING_SEED_FC_HEIGHT=720
CRAWLING_SEED_FC_FRAMERATE=10
# Each render worker may load the dense dconn; keep this low unless RAM is ample.
CRAWLING_SEED_FC_NUM_CPUS=1
CRAWLING_SEED_FC_TARGET_SIZE_MB=10
CRAWLING_SEED_FC_FORCE=0
CRAWLING_SEED_FC_FORCE_DCONN=0
CRAWLING_SEED_FC_SKIP_MOVIE=0
CRAWLING_SEED_FC_KEEP_DCONN=0

# -----------------------------------------------------------------------------
# 3g) NSI
# -----------------------------------------------------------------------------
NSI_USE_EXTERNAL_CLI=1
NSI_PYTHON="$CONCAT_PYTHON"
NSI_INPUT_TAG=""                 # empty => derive from CONCAT_INPUT_TAG
NSI_CONCAT_OUT_SUBDIR="$CONCAT_OUT_SUBDIR"
NSI_FD_THRESHOLD="$CONCAT_FD_THRESHOLD"
NSI_USABILITY_MODEL=1
NSI_RELIABILITY_MODEL=0

# -----------------------------------------------------------------------------
# 3h) PFM
# -----------------------------------------------------------------------------
PFM_STRATEGY="ridge_fusion"      # ridge_fusion | infomap
PFM_PYTHON="$PIPELINE_PYTHON"
PFM_INPUT_CIFTI=""               # empty => derive from concat outputs below
PFM_INPUT_TAG=""                 # empty => derive from CONCAT_INPUT_TAG
PFM_CONCAT_OUT_SUBDIR="$CONCAT_OUT_SUBDIR"
PFM_FD_THRESHOLD="$CONCAT_FD_THRESHOLD"
PFM_DISTANCE_MATRIX=""           # empty => <SubjectDir>/anat/T1w/fsaverage_LR32k/DistanceMatrix.npy
PFM_DISTANCE_BUILD_IF_MISSING=1
PFM_DISTANCE_CORTEX_MODE="hybrid"      # geodesic | hybrid | euclidean
PFM_DISTANCE_EUCLIDEAN_OVERRIDE_MM=5 # hybrid only: use 3D Euclidean distance for same-hemi cortical pairs this close
PFM_OUTDIR=""                    # empty => <SubjectDir>/func/<FUNC_DIRNAME>/PFM
PFM_PREP_DIR=""                  # empty => <PFM_OUTDIR>/prep

# Optional post-PFM manual update pass. This scans per-density review tables and
# updates labeled density and consensus outputs without rerunning Infomap.
PFM_INFOMAP_UPDATE_ENABLE=0
PFM_INFOMAP_UPDATE_TABLE_GLOB="GraphDensity_*/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv"
PFM_INFOMAP_UPDATE_OUTFILE="InfomapNetworkLabels_ManualAdjusted"

# Optional homogeneity/spin-null diagnostics. A rotation-index CIFTI should use
# 1-based source indices in cortex rows, with 0 for unmapped vertices.
PFM_HOMOGENEITY_TEST_ENABLE=0
PFM_HOMOGENEITY_ROTATIONS_CIFTI="$MEDIR/res0urces/Rotated_inds.dtseries.nii"
PFM_HOMOGENEITY_N_ROTATIONS=100
PFM_HOMOGENEITY_MIN_COMMUNITY_SIZE=5
PFM_HOMOGENEITY_MAX_MEMBERS_PER_COMMUNITY=1000
PFM_HOMOGENEITY_ALPHA=0.05

# =============================================================================
# 4) Advanced / Rarely Changed Knobs
# =============================================================================

# -----------------------------------------------------------------------------
# 4a) Module Entrypoints
# -----------------------------------------------------------------------------
VALIDATE_MODULE="$MEDIR/modules/mefmri_validate_inputs.sh"
ANAT_HCP_MODULE="$MEDIR/modules/mefmri_anat_hcp.sh"
ANAT_CHARM_MODULE="$MEDIR/modules/mefmri_anat_charm.sh"
FUNC_FIELDMAPS_MODULE="$MEDIR/modules/mefmri_func_fieldmaps.sh"
FUNC_MEDIC_MODULE="$MEDIR/modules/mefmri_func_medic.sh"
FUNC_COREG_MODULE="$MEDIR/modules/mefmri_func_coreg.sh"
FUNC_HEADMOTION_MODULE="$MEDIR/modules/mefmri_func_headmotion.sh"
FUNC_MEICA_MODULE="$MEDIR/modules/mefmri_func_meica.sh"
FUNC_SINGLEECHO_MODULE="$MEDIR/modules/mefmri_func_singleecho.sh"
FUNC_ACOMPCOR_MODULE="$MEDIR/modules/mefmri_func_acompcor.sh"
FUNC_AROMA_MODULE="$MEDIR/modules/mefmri_func_aroma.sh"
FUNC_MGTR_MODULE="$MEDIR/modules/mefmri_func_mgtr.sh"
FUNC_VOL2SURF_MODULE="$MEDIR/modules/mefmri_func_vol2surf.sh"
FUNC_CONCAT_MODULE="$MEDIR/modules/mefmri_func_concat.sh"
FUNC_NSI_MODULE="$MEDIR/modules/mefmri_func_nsi.sh"
FUNC_PFM_MODULE="$MEDIR/modules/mefmri_func_pfm.sh"
FUNC_PFM_INFOMAP_UPDATE_MODULE="$MEDIR/modules/mefmri_func_pfm_infomap_update.sh"
FUNC_CRAWLING_SEED_FC_MODULE="$MEDIR/modules/mefmri_func_crawling_seed_fc.sh"

# -----------------------------------------------------------------------------
# 4b) Vol2Surf GoodVoxels / Masking Details
# -----------------------------------------------------------------------------
VOL2SURF_USE_GOOD_VOXELS_MASK=0
VOL2SURF_GOOD_VOXELS_FACTOR="0.5"
VOL2SURF_GOOD_VOXELS_SIGMA_MM="5"
VOL2SURF_GOOD_VOXELS_KEEP_INTERMEDIATES=1

# -----------------------------------------------------------------------------
# 4c) MEDIC / warpkit Detailed Options
# -----------------------------------------------------------------------------
# Experimental options used only when DISTORTION_CORRECTION_MODE="medic".
MEDIC_OVERWRITE=0                 # 0 skips existing MEDIC outputs, 1 reruns wk-medic
MEDIC_NOISE_FRAMES=0              # trailing frames to trim before MEDIC, passed to wk-medic --noiseframes
MEDIC_DEBUG=0                     # 1 passes --debug to wk-medic
MEDIC_WRAP_LIMIT=0                # 1 passes --wrap-limit to wk-medic
MEDIC_COMPUTE_JACOBIAN=1          # 1 also writes per-frame Jacobian maps

# -----------------------------------------------------------------------------
# 4d) Tedana Detailed Options
# -----------------------------------------------------------------------------
TEDANA_FITTYPE="curvefit"        # curvefit|loglin
TEDANA_ICA_METHOD="fastica"      # fastica|robustica
TEDANA_N_ROBUST_RUNS=10
TEDANA_SEED=42
TEDANA_THREADS=4
TEDANA_MASKTYPE="none"           # none|dropout|decay
TEDANA_CONVENTION="orig"         # orig|bids
TEDANA_OVERWRITE=1
TEDANA_LOWMEM=0
TEDANA_USE_EXTERNAL_MIX=0
TEDANA_EXTERNAL_MIX_BASENAME=""

# -----------------------------------------------------------------------------
# 4e) MEICA / Reclassify Detailed Thresholds
# -----------------------------------------------------------------------------
MEICA_PRIORS_MAT="$PFM_NETWORK_PRIORS_CORTICAL_MAT"
MEICA_BETAS_CIFTI=""
MEICA_RHO_RESCUE=0.30
MEICA_RHO_REJECT=0.10

MEICA_NSI_RESCUE_THRESHOLD=0.20
MEICA_NSI_RESCUE_QUANTILE=0.10
MEICA_NSI_KILL_MODE="adaptive"   # adaptive|fixed
MEICA_NSI_KILL_THRESHOLD=0.05    # used only when kill mode=fixed
MEICA_NSI_KILL_MIN=0.04
MEICA_NSI_KILL_MAX=0.10
MEICA_NSI_KILL_INTERCEPT=0.14
MEICA_NSI_KILL_SLOPE=0.25
MEICA_NSI_GUARDRAIL_KAPPA_RHO=1
MEICA_SUBCORT_RATIO_THRESH=5.0
MEICA_KILL_PRIORITY_ENABLE=0
MEICA_KILL_PRIORITY_W_LOGRATIO=0.50
MEICA_KILL_PRIORITY_W_NSI=0.30
MEICA_KILL_PRIORITY_W_VAR=0.20
MEICA_KILL_VAR_FLOOR_QUANTILE=0.60
MEICA_KILL_CUMVAR_CAP=0.95

# -----------------------------------------------------------------------------
# 4f) Concat Detailed Options
# -----------------------------------------------------------------------------
CONCAT_DEMEAN_RUNS=1
CONCAT_VAR_NORM_RUNS=0
CONCAT_VAR_NORM_EPS=1e-8
CONCATENATE_RUNS=1
CONCAT_SAVE_FD_TXT=1
CONCAT_SAVE_SCANIDX_TXT=1

# -----------------------------------------------------------------------------
# 4g) NSI Detailed Options
# -----------------------------------------------------------------------------
NSI_RELIABILITY_NSI_T=10
NSI_RELIABILITY_QUERY_T=60
NSI_EXTERNAL_ROOT="$MEDIR/lib/pfm-nsi"
NSI_EXTERNAL_ENTRY="pfm_nsi.cli"
NSI_EXTERNAL_OUT_SUBDIR=""              # empty => write directly to func/<FUNC_DIRNAME>/qa/NSI
NSI_EXTERNAL_PREFIX="pfm_nsi"
NSI_EXTERNAL_USABILITY="$NSI_USABILITY_MODEL"
NSI_EXTERNAL_RELIABILITY="$NSI_RELIABILITY_MODEL"
NSI_EXTERNAL_NSI_T="$NSI_RELIABILITY_NSI_T"
NSI_EXTERNAL_QUERY_T="$NSI_RELIABILITY_QUERY_T"
NSI_EXTERNAL_THRESHOLDS="0.6,0.7,0.8"
NSI_EXTERNAL_STRUCTURES=""              # empty => pfm-nsi default bilateral structures
NSI_EXTERNAL_MORANS=0
NSI_EXTERNAL_SLOPE=0
NSI_EXTERNAL_RIDGE_LAMBDAS="10"
NSI_EXTERNAL_SPARSE_FRAC=""
NSI_EXTERNAL_THREADS=4
NSI_EXTERNAL_FULLMEM=0
NSI_EXTERNAL_DTYPE="float32"
NSI_EXTERNAL_BLOCK_SIZE=2048
NSI_EXTERNAL_KEEP_ALLRHO=1
NSI_EXTERNAL_KEEP_BETAS=1
NSI_EXTERNAL_KEEP_FC_MAP=0

# -----------------------------------------------------------------------------
# 4h) PFM Shared Resources
# -----------------------------------------------------------------------------
PFM_DISTANCE_VARIANT_CHUNK_ROWS=128
PFM_RESOURCES_ROOT="$MEDIR/res0urces"

# -----------------------------------------------------------------------------
# 4i) PFM Ridge-Fusion
# -----------------------------------------------------------------------------
PFM_RF_OUTFILE="RidgeFusion_VTX"
PFM_RF_FC_WEIGHT=1.0              # weight on subject functional connectivity evidence
PFM_RF_FC_DEMEAN=0             # set 1 when ridge_fusion input has not had GSR/MGTR; demeans each FC fingerprint
PFM_RF_SPATIAL_WEIGHT=0.2         # weight on spatial/network priors
PFM_RF_LAMBDA=10                  # ridge penalty; larger values shrink estimates more strongly toward priors
PFM_RF_LOCAL_EXCLUSION_MM=30      # excludes nearby vertices from FC fingerprints to avoid local smoothing bias
PFM_RF_SUBCORT_REGRESS_ENABLE=1   # regress nearby cortical signal from subcortex before PFM estimation
PFM_RF_SUBCORT_REGRESS_DISTANCE_MM=20
PFM_RF_BRAIN_STRUCTURES_CSV="CORTEX_LEFT,CEREBELLUM_LEFT,ACCUMBENS_LEFT,CAUDATE_LEFT,PUTAMEN_LEFT,THALAMUS_LEFT,HIPPOCAMPUS_LEFT,AMYGDALA_LEFT,CORTEX_RIGHT,CEREBELLUM_RIGHT,ACCUMBENS_RIGHT,CAUDATE_RIGHT,PUTAMEN_RIGHT,THALAMUS_RIGHT,HIPPOCAMPUS_RIGHT,AMYGDALA_RIGHT"
PFM_RF_SMOOTHING_KERNEL=1.7       # mm; set 0 to disable CIFTI smoothing before PFM estimation

# Optional areal parcellation of the network-level PFM output.
PFM_AREAL_ENABLE=1
PFM_AREAL_OUTFILE=""              # empty => auto: <network-label-output>+ArealParcellation
PFM_AREAL_MIN_SIZE=30             # minimum parcel size in grayordinates
PFM_PRIORS_MAT="$PFM_NETWORK_PRIORS_CORTICAL_MAT"
PFM_SUBCORT_PRIORS_NII="$PFM_NETWORK_PRIORS_SUBCORTICAL_NII"
PFM_NEIGHBORS_MAT="$PFM_RESOURCES_ROOT/Cifti_surf_neighbors_LR_normalwall.mat"

# Experimental multi-network Ridge-Fusion extension (requires PFM_AREAL_ENABLE=1).
# Disabled by default; emits null-selected whole-parcel and localized maps.
# Preliminary behavior is encouraging; keep opt-in pending broader validation.
PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE=0
PFM_RF_MULTINETWORK_STRIPE_PRESET="balanced" # strict | balanced | loose | custom
PFM_RF_MULTINETWORK_OUTDIR_NAME="ExperimentalMultiNetwork"
PFM_RF_MULTINETWORK_REFERENCE_COUNT=1500
PFM_RF_MULTINETWORK_REFERENCE_SEED=12345
PFM_RF_MULTINETWORK_CV_LOCAL_EXCLUSION_MM=30
PFM_RF_MULTINETWORK_SPLIT_REPEATS=8
PFM_RF_MULTINETWORK_SPLIT_BLOCK_LENGTH=25
PFM_RF_MULTINETWORK_SPLIT_SEED=12345
# Presets override consistency, correction, weights, and local thresholds below;
# use STRIPE_PRESET="custom" to control those individually. Seeds/repeats remain active.
PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION=0.8
PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION=0.8
PFM_RF_MULTINETWORK_SELECTION_CORRECTION="fdr" # fdr (descriptive default) | fwer (extreme-confidence tier)
PFM_RF_MULTINETWORK_FDR_ALPHA=0.05
PFM_RF_MULTINETWORK_PARCEL_NULL_ITERATIONS=2000
PFM_RF_MULTINETWORK_PARCEL_NULL_SEED=12345
PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT=0.20
PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT=0.15
PFM_RF_MULTINETWORK_NULL_ALPHA=0.05
PFM_RF_MULTINETWORK_MAX_NETWORKS=3
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN=0.10
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN=0.40
PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES=20
PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS=5

# -----------------------------------------------------------------------------
# 4j) PFM Infomap
# -----------------------------------------------------------------------------
PFM_INFOMAP_DISTANCE_MATRIX=""          # empty => PFM_DISTANCE_MATRIX (.npy; auto-built when missing)
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001" # graph densities, from dense to sparse
PFM_INFOMAP_NUM_REPS_EXPR="5,10,20,30,50,75,100"                              # one repeat count per density
PFM_INFOMAP_MIN_DISTANCE=10       # mm; suppresses local edges before community detection
PFM_INFOMAP_BAD_VERTS_CSV=""      # optional comma-separated 0-based grayordinate indices to exclude
PFM_INFOMAP_STRUCTURES_CSV="CORTEX_LEFT,CEREBELLUM_LEFT,ACCUMBENS_LEFT,CAUDATE_LEFT,PUTAMEN_LEFT,THALAMUS_LEFT,HIPPOCAMPUS_LEFT,CORTEX_RIGHT,CEREBELLUM_RIGHT,ACCUMBENS_RIGHT,CAUDATE_RIGHT,PUTAMEN_RIGHT,THALAMUS_RIGHT,HIPPOCAMPUS_RIGHT"
PFM_INFOMAP_NUM_CORES=1
PFM_INFOMAP_BINARY=""                   # set explicit path if infomap is not already on PATH
PFM_INFOMAP_DRY_RUN=0                   # 1=validate arguments/wiring without running computation
PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1    # 1 labels communities with canonical prior network IDs
PFM_INFOMAP_LABEL_OUTFILE="InfomapNetworkLabels"
PFM_INFOMAP_LABEL_FC_WEIGHT=1.0         # community-to-prior functional similarity weight
PFM_INFOMAP_LABEL_SPATIAL_WEIGHT=1.0    # community-to-prior spatial overlap weight
PFM_INFOMAP_LABEL_CONFIDENCE_THRESHOLD=0.15 # below this margin, labels are written but flagged for review
PFM_INFOMAP_LABEL_MIN_FC_SIMILARITY=0.33
PFM_INFOMAP_LABEL_MIN_COMMUNITY_SIZE=10
PFM_INFOMAP_LABEL_UNASSIGNED_VALUE=21   # label ID for communities too small or too weak to assign
PFM_INFOMAP_LABEL_STRICT_THRESHOLDING=0  # 1 converts low-confidence communities to unassigned
PFM_INFOMAP_LABEL_DENSITY_INDEX=-1       # -1 labels all density columns; positive values are 1-based
PFM_INFOMAP_MANUAL_LABEL_APPLY_ENABLE=0
PFM_INFOMAP_MANUAL_LABEL_TABLE=""        # empty => <PFM_OUTDIR>/<PFM_INFOMAP_LABEL_OUTFILE>_ManualCorrections.csv
PFM_INFOMAP_MANUAL_LABEL_OUTFILE="InfomapNetworkLabels_ManualAdjusted"

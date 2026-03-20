#!/usr/bin/env bash
# Config for RevisedMe-fMRIPipeline/mefmri_pipeline.sh
#
# Wrapper call:
#   mefmri_pipeline.sh <SubjectDir> [ConfigFile]
#
# =============================================================================
# 1) Required Paths and Environment Settings
# =============================================================================
# TODO(user): review and set machine-specific values in this section before
# first run on a new system. At minimum, confirm:
# - CHARM_BIN (if charm is not on PATH)
# - TEDANA_ENV / TEDANA_ACTIVATE_MODE
# - NSI_EXTERNAL_ROOT (only if NSI_USE_EXTERNAL_CLI=1)
# - FS_LICENSE / FS_LICENSE_FILE in your shell environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EnvironmentScript="$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh"

# TODO(user): set CHARM_BIN to your local SimNIBS CHARM binary path if "charm"
# is not already on PATH. Leave empty to use PATH lookup in modules.
CHARM_BIN=""

# PFM network priors:
# - cortical prior weights / labels in MATLAB format
# - subcortical prior volume in NIfTI format
PFM_NETWORK_PRIORS_CORTICAL_MAT="$MEDIR/res0urces/Priors.mat"
PFM_NETWORK_PRIORS_SUBCORTICAL_NII="$MEDIR/res0urces/SubcorticalPriors.nii.gz"

# Threading model 
THREADS_DEFAULT=8
THREADS_ANAT_HCP="$THREADS_DEFAULT"
THREADS_FIELDMAPS="$THREADS_DEFAULT"
THREADS_COREG="$THREADS_DEFAULT"
THREADS_HEADMOTION="$THREADS_DEFAULT"
THREADS_MEICA=2

# Module entrypoint mapping (rename-safe)
VALIDATE_MODULE="$MEDIR/modules/mefmri_validate_inputs.sh"
ANAT_HCP_MODULE="$MEDIR/modules/mefmri_anat_hcp.sh"
ANAT_CHARM_MODULE="$MEDIR/modules/mefmri_anat_charm.sh"
FUNC_FIELDMAPS_MODULE="$MEDIR/modules/mefmri_func_fieldmaps.sh"
FUNC_COREG_MODULE="$MEDIR/modules/mefmri_func_coreg.sh"
FUNC_HEADMOTION_MODULE="$MEDIR/modules/mefmri_func_headmotion.sh"
FUNC_MEICA_MODULE="$MEDIR/modules/mefmri_func_meica.sh"
FUNC_MGTR_MODULE="$MEDIR/modules/mefmri_func_mgtr.sh"
FUNC_VOL2SURF_MODULE="$MEDIR/modules/mefmri_func_vol2surf.sh"
FUNC_CONCAT_MODULE="$MEDIR/modules/mefmri_func_concat.sh"
FUNC_NSI_MODULE="$MEDIR/modules/mefmri_func_nsi.sh"
FUNC_PFM_MODULE="$MEDIR/modules/mefmri_func_pfm.sh"

# =============================================================================
# 2) Common User-Facing Pipeline Options
# =============================================================================
# Resume controls
START_SESSION=1
START_FROM_MODULE="validate"   # validate|anat_hcp|anat_charm|fieldmaps|coreg|headmotion|meica|mgtr|vol2surf|concat|nsi|pfm
STOP_AFTER_MODULE="pfm"   # "" to run full chain

# Shared routing and naming
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"
# Internal transform-space folder used by coreg outputs (kept as "rest" for compatibility).
FUNC_XFMS_DIRNAME="rest"
DOF=6
AtlasTemplate="$MEDIR/res0urces/MNI152_T1_2mm.nii.gz"
AtlasSpace="T1w"               # T1w|MNINonlinear
APPLY_N4_BIAS=0                 # 0|1 ; If PreScan Normalize is turned off at acquistion, you cna probably leave this off.

# Pipeline-level output/provenance
RUN_CONFIG_SNAPSHOT=1           # 1: write effective run metadata snapshot into subject func/qa

# =============================================================================
# 3) Masking and Surface-Mapping Settings
# =============================================================================
# CHARM anatomical mask mode used to regenerate T1/T2 brain images.
# - charm: CHARM labeling-derived mask (<100), then dilated [default]
# - hcp: preserve pre-CHARM HCP-style whole-brain mask
CHARM_BRAIN_MASK_MODE="charm"  # charm|hcp
CHARM_BRAIN_MASK_DILATE_ITERS=1 # integer >= 0 (used in charm mode)

# Cortical ribbon generation in anat module
CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS=1  # 1 removes HC/Amy/CSF, 0 keeps FS ribbon-only labels
CHARM_WRITE_CORTICAL_RIBBON=1           # 1 writes anat/T1w/CorticalRibbon.nii.gz

# MGTR cortical ribbon source
# - xfms: func/xfms/rest/CorticalRibbon_*_func_mask.nii.gz
# - legacy_rois: func/rois/CorticalRibbon.nii.gz built from FreeSurfer ribbon mgz
MGTR_RIBBON_SOURCE="xfms"      # xfms|legacy_rois

# Vol2Surf inputs and mapping behavior
VOL2SURF_INPUTS="OCME,OCME+MEICA,OCME+MEICA+MGTR"
VOL2SURF_USE_CORTICAL_RIBBON_MASK=1    # Default modern behavior
VOL2SURF_CIFTI_STAMP=""               # Optional non-canonical suffix for experimental comparisons

# =============================================================================
# 4) Tedana and Denoising Settings
# =============================================================================
# Modern defaults:
# - tedana 0.26.00 path (via environment)
# - Kundu PCA selection
TEDANA_ENV="mefmri_env" 			  # name of the conda environment used for tedana
TEDANA_ACTIVATE_MODE="conda_activate" # conda_activate|conda_run|direct
TEDANA_COMPAT_MODE="modern"           # recommended default
MEPCA="350"                           # default kundu

MaxIterations=500
MaxRestarts=5
TEDANA_FITTYPE="curvefit"             # curvefit|loglin
TEDANA_ICA_METHOD="fastica"           # fastica|robustica
TEDANA_N_ROBUST_RUNS=10
TEDANA_SEED=42
TEDANA_THREADS=4
TEDANA_MASKTYPE="none"                # none|dropout|decay
TEDANA_CONVENTION="orig"              # orig|bids
TEDANA_OVERWRITE=1
TEDANA_LOWMEM=0
TEDANA_USE_EXTERNAL_MIX=0
TEDANA_EXTERNAL_MIX_BASENAME=""
MEICA_SKIP_TEDANA_IF_EXISTS=0

# =============================================================================
# 5) Reclassification Settings
# =============================================================================
MEICA_PARALLEL_JOBS=4
MEICA_RECLASSIFY_ENABLE=1
MEICA_TEDANA_SUBDIR="Tedana"
MEICA_RECLASS_SUBDIR="Reclassify"

# Default classifier is modern NSI-based mode.
# Optional fallback: legacy_template_rho (see section 9 for legacy-focused presets)
MEICA_CLASSIFIER_MODE="nsi"           # nsi|legacy_template_rho|none
MEICA_PRIORS_MAT="$PFM_NETWORK_PRIORS_CORTICAL_MAT"
MEICA_BETAS_CIFTI=""
MEICA_RHO_RESCUE=0.30
MEICA_RHO_REJECT=0.10

# Adaptive NSI kill defaults
MEICA_NSI_RESCUE_THRESHOLD=0.20
MEICA_NSI_RESCUE_QUANTILE=0.10
MEICA_NSI_KILL_MODE="adaptive"        # adaptive|fixed
MEICA_NSI_KILL_THRESHOLD=0.05         # used only when kill mode=fixed
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

# =============================================================================
# 6) Output, Reporting, and Overwrite Behavior
# =============================================================================
# Reclassify HTML reports are regenerated by default.
MEICA_RECLASS_NO_REPORTS=0             # 0: write reports (default), 1: suppress report generation
MEICA_QC_CIFTI_ENABLE=1
MEICA_QC_CIFTI_TAGS="betas_OC,t2sv,s0v"
MEICA_ORIG_ALIAS_ENABLE=1             # 1: create legacy/orig filename aliases from modern tedana desc-* outputs

# HCP anat module behavior
HCP_ANAT_CLEAN_START=1                 # 1 removes prior anat outputs before rerun
CHARM_REUSE_EXISTING_M2M=0             # 1 reuses existing anat/m2m_<Subject> outputs if present

# =============================================================================
# 7) Atlas / Template / Reference Inputs
# =============================================================================
# Fieldmap/Coreg metadata options
FM_PE_MODE="infer"                    # infer|config
FM_AP_PE_DIR=""                       # e.g. j-
FM_PA_PE_DIR=""                       # e.g. j
EPIREG_PEDIR=""                       # empty => infer from PE.txt
SCAN_SPECIFIC_FM=1
HEADMOTION_KEEP_MCF=0

# =============================================================================
# 8) Concat Module Settings
# =============================================================================
# Run-level concatenation + FD censor prep.
CONCAT_ENABLE=1
CONCAT_PYTHON="python3"
CONCAT_INPUT_TAG="OCME+MEICA+MGTR"
CONCAT_OUT_SUBDIR="ConcatenatedCiftis"
CONCAT_DEMEAN_RUNS=1                    # demean each run before concatenation
CONCAT_VAR_NORM_RUNS=0                  # optional run-wise variance normalization after demeaning
CONCAT_VAR_NORM_EPS=1e-8                # floor to avoid divide-by-zero during variance normalization
CONCATENATE_RUNS=1
CONCAT_CENSOR_BY_FD=1
CONCAT_FD_THRESHOLD=0.3
CONCAT_SAVE_FD_TXT=1
CONCAT_SAVE_SCANIDX_TXT=1

# =============================================================================
# 9) NSI Module Settings
# =============================================================================
# Can run independently after concat.
NSI_ENABLE=1
NSI_USE_EXTERNAL_CLI=1
NSI_PYTHON="$CONCAT_PYTHON"
NSI_INPUT_TAG="$CONCAT_INPUT_TAG"
NSI_CONCAT_OUT_SUBDIR="$CONCAT_OUT_SUBDIR"
NSI_FD_THRESHOLD="$CONCAT_FD_THRESHOLD"
NSI_USABILITY_MODEL=1
NSI_RELIABILITY_MODEL=0
NSI_RELIABILITY_NSI_T=10
NSI_RELIABILITY_QUERY_T=60
# Bundled local NSI fork.
NSI_EXTERNAL_ROOT="$MEDIR/lib/pfm-nsi"
NSI_EXTERNAL_ENTRY="pfm_nsi.cli"
NSI_EXTERNAL_OUT_SUBDIR=""              # empty => write directly to func/qa/NSI
NSI_EXTERNAL_PREFIX="pfm_nsi"
NSI_EXTERNAL_USABILITY="$NSI_USABILITY_MODEL"     # compatibility alias
NSI_EXTERNAL_RELIABILITY="$NSI_RELIABILITY_MODEL" # compatibility alias
NSI_EXTERNAL_NSI_T="$NSI_RELIABILITY_NSI_T"       # compatibility alias
NSI_EXTERNAL_QUERY_T="$NSI_RELIABILITY_QUERY_T"   # compatibility alias
NSI_EXTERNAL_THRESHOLDS="0.6,0.7,0.8"
NSI_EXTERNAL_STRUCTURES=""               # empty => pfm-nsi default bilateral structures; otherwise comma-separated CIFTI names
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

# =============================================================================
# 10) PFM Module Settings
# =============================================================================
PFM_ENABLE=1
PFM_STRATEGY="ridge_fusion"             # ridge_fusion | infomap
PFM_PYTHON="python3"

# Inputs from prior stages
PFM_INPUT_CIFTI=""                      # empty => derive from concat outputs below
PFM_INPUT_TAG="$CONCAT_INPUT_TAG"
PFM_CONCAT_OUT_SUBDIR="$CONCAT_OUT_SUBDIR"
PFM_FD_THRESHOLD="$CONCAT_FD_THRESHOLD"
PFM_DISTANCE_MATRIX=""                  # empty => <SubjectDir>/anat/T1w/fsaverage_LR32k/DistanceMatrix.npy
PFM_DISTANCE_VARIANT_CHUNK_ROWS=128
PFM_DISTANCE_BUILD_IF_MISSING=1
PFM_OUTDIR=""                           # empty => <SubjectDir>/func/<FUNC_DIRNAME>/PFM

# External PFM resource roots
PFM_RESOURCES_ROOT="$MEDIR/res0urces"
PFM_INFOMAP_WRAPPER="$PFM_RESOURCES_ROOT/PFM-InfoMap-Tmp/pfm_wrapper.m"

# Ridge-fusion 
# strategy settings
PFM_RF_ENABLE=1
PFM_RF_OUTFILE="RidgeFusion_VTX"
PFM_RF_FC_WEIGHT=1.0
PFM_RF_SPATIAL_WEIGHT=0.1
PFM_RF_LAMBDA=10
PFM_RF_LOCAL_EXCLUSION_MM=10
PFM_RF_SUBCORT_REGRESS_ENABLE=1
PFM_RF_SUBCORT_REGRESS_DISTANCE_MM=20
PFM_RF_BRAIN_STRUCTURES_CSV="CORTEX_LEFT,CEREBELLUM_LEFT,ACCUMBENS_LEFT,CAUDATE_LEFT,PUTAMEN_LEFT,THALAMUS_LEFT,HIPPOCAMPUS_LEFT,AMYGDALA_LEFT,CORTEX_RIGHT,CEREBELLUM_RIGHT,ACCUMBENS_RIGHT,CAUDATE_RIGHT,PUTAMEN_RIGHT,THALAMUS_RIGHT,HIPPOCAMPUS_RIGHT,AMYGDALA_RIGHT"
PFM_RF_SMOOTHING_KERNEL=1.7             # mm; <=0 disables smoothing step

# Reserved knobs for downstream areal parcellation and prior controls
PFM_AREAL_ENABLE=1
PFM_AREAL_OUTFILE="RidgeFusion_VTX+ArealParcellation"
PFM_AREAL_MIN_SIZE=10
PFM_AP_METHOD="kmeans"
PFM_AP_KMAX_CAP=8
PFM_AP_KMIN_PER_PATCH=1
PFM_AP_DIFFUSE_ITERS=6
PFM_AP_LAMBDA_SPATIAL=8
PFM_PRIORS_MAT="$PFM_NETWORK_PRIORS_CORTICAL_MAT"
PFM_SUBCORT_PRIORS_NII="$PFM_NETWORK_PRIORS_SUBCORTICAL_NII"
PFM_NEIGHBORS_MAT="$PFM_RESOURCES_ROOT/Cifti_surf_neighbors_LR_normalwall.mat"

# Infomap strategy settings 
PFM_INFOMAP_DISTANCE_MATRIX=""          # empty => <SubjectDir>/anat/T1w/fsaverage_LR32k/DistanceMatrix.mat
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="0.001:0.001:0.009"
PFM_INFOMAP_NUM_REPS_EXPR="20"
PFM_INFOMAP_MIN_DISTANCE=20
PFM_INFOMAP_BAD_VERTS_EXPR="[]"
PFM_INFOMAP_STRUCTURES_EXPR="[]"
PFM_INFOMAP_BAD_VERTS_CSV=""
PFM_INFOMAP_STRUCTURES_CSV=""
PFM_INFOMAP_NUM_CORES=1
PFM_INFOMAP_BINARY=""                   # set explicit path if infomap is not already on PATH
PFM_INFOMAP_USE_MATLAB=0                # 0=Python default; 1=MATLAB pfm_wrapper fallback
PFM_INFOMAP_DRY_RUN=0                   # 1=validate arguments/wiring without running computation
PFM_INFOMAP_NETWORK_MAPPING_ENABLE=0    # keep 0 until mapping strategy is finalized

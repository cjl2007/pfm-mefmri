#!/usr/bin/env bash
# HCP anatomical preprocessing wrapper (PreFreeSurfer + FreeSurfer + PostFreeSurfer).
set -euo pipefail

StudyFolder="${1:?Need StudyFolder}"
Subject="${2:?Need Subject}"
export NSLOTS="${3:?Need NTHREADS}"

if [ "${StudyFolder: -1}" = "/" ]; then
  StudyFolder="${StudyFolder%?}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="${MEDIR:-$SCRIPT_DIR}"
EnvironmentScript="${EnvironmentScript:-$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh}"
PRINTCOM="${PRINTCOM:-}"

# Enable this only when you explicitly want to clear previous anat outputs.
HCP_ANAT_CLEAN_START="${HCP_ANAT_CLEAN_START:-0}"

# Distortion/fmap defaults.
AvgrdcSTRING="${AvgrdcSTRING:-NONE}"
MagnitudeInputName="${MagnitudeInputName:-NONE}"
PhaseInputName="${PhaseInputName:-NONE}"
TE="${TE:-NONE}"
SpinEchoPhaseEncodeNegative="${SpinEchoPhaseEncodeNegative:-NONE}"
SpinEchoPhaseEncodePositive="${SpinEchoPhaseEncodePositive:-NONE}"
SEEchoSpacing="${SEEchoSpacing:-NONE}"
SEUnwarpDir="${SEUnwarpDir:-NONE}"
TopupConfig="${TopupConfig:-NONE}"
GEB0InputName="${GEB0InputName:-NONE}"

if [ ! -f "$EnvironmentScript" ]; then
  echo "ERROR: missing HCP environment script: $EnvironmentScript" >&2
  exit 2
fi
source "$EnvironmentScript"

T1wTemplate="${HCPPIPEDIR_Templates}/MNI152_T1_0.8mm.nii.gz"
T1wTemplateBrain="${HCPPIPEDIR_Templates}/MNI152_T1_0.8mm_brain.nii.gz"
T1wTemplate2mm="${HCPPIPEDIR_Templates}/MNI152_T1_2mm.nii.gz"
T2wTemplate="${HCPPIPEDIR_Templates}/MNI152_T2_0.8mm.nii.gz"
T2wTemplateBrain="${HCPPIPEDIR_Templates}/MNI152_T2_0.8mm_brain.nii.gz"
T2wTemplate2mm="${HCPPIPEDIR_Templates}/MNI152_T2_2mm.nii.gz"
TemplateMask="${HCPPIPEDIR_Templates}/MNI152_T1_0.8mm_brain_mask.nii.gz"
Template2mmMask="${HCPPIPEDIR_Templates}/MNI152_T1_2mm_brain_mask_dil.nii.gz"

T1wSampleSpacing="${T1wSampleSpacing:-NONE}"
T2wSampleSpacing="${T2wSampleSpacing:-NONE}"
UnwarpDir="${UnwarpDir:-z}"
BrainSize="${BrainSize:-170}"
FNIRTConfig="${FNIRTConfig:-${HCPPIPEDIR_Config}/T1_2_MNI152_2mm.cnf}"
GradientDistortionCoeffs="${GradientDistortionCoeffs:-NONE}"

ANAT_BASE="${StudyFolder}/${Subject}/anat"
UNPROC_T1_DIR="${ANAT_BASE}/unprocessed/T1w"
UNPROC_T2_DIR="${ANAT_BASE}/unprocessed/T2w"
QA_DIR="${StudyFolder}/${Subject}/qa"

if [ ! -d "$UNPROC_T1_DIR" ]; then
  echo "ERROR: missing T1w input directory: $UNPROC_T1_DIR" >&2
  exit 2
fi

if [[ "$HCP_ANAT_CLEAN_START" == "1" ]]; then
  rm -rf "${StudyFolder}/${Subject}/T1w" "${StudyFolder}/${Subject}/T2w" \
    "${StudyFolder}/${Subject}/MNINonLinear" "${StudyFolder}/${Subject}/qa" \
    "${ANAT_BASE}/T1w" "${ANAT_BASE}/T2w" "${ANAT_BASE}/MNINonLinear" "${ANAT_BASE}/qa"
fi

T1wInputImages=""
for i in "${UNPROC_T1_DIR}"/T1w*.nii.gz; do
  [ -e "$i" ] || continue
  T1wInputImages="${T1wInputImages}${i}@"
done

if [ -z "$T1wInputImages" ]; then
  echo "ERROR: no T1w inputs found in ${UNPROC_T1_DIR}" >&2
  exit 2
fi

T2wInputImages=""
for i in "${UNPROC_T2_DIR}"/T2w*.nii.gz; do
  [ -e "$i" ] || continue
  T2wInputImages="${T2wInputImages}${i}@"
done

if [ -z "$T2wInputImages" ]; then
  T2wInputImages="NONE"
  ProcessingMode="LegacyStyleData"
else
  ProcessingMode="HCPStyleData"
fi

mkdir -p "$QA_DIR"

echo "Running PreFreeSurferPipeline for ${Subject}"
"${HCPPIPEDIR}/PreFreeSurfer/PreFreeSurferPipeline.sh" \
  --path="$StudyFolder" \
  --subject="$Subject" \
  --t1="$T1wInputImages" \
  --t2="$T2wInputImages" \
  --t1template="$T1wTemplate" \
  --t1templatebrain="$T1wTemplateBrain" \
  --t1template2mm="$T1wTemplate2mm" \
  --t2template="$T2wTemplate" \
  --t2templatebrain="$T2wTemplateBrain" \
  --t2template2mm="$T2wTemplate2mm" \
  --templatemask="$TemplateMask" \
  --template2mmmask="$Template2mmMask" \
  --brainsize="$BrainSize" \
  --fnirtconfig="$FNIRTConfig" \
  --fmapmag="$MagnitudeInputName" \
  --fmapphase="$PhaseInputName" \
  --fmapgeneralelectric="$GEB0InputName" \
  --echodiff="$TE" \
  --SEPhaseNeg="$SpinEchoPhaseEncodeNegative" \
  --SEPhasePos="$SpinEchoPhaseEncodePositive" \
  --seechospacing="$SEEchoSpacing" \
  --seunwarpdir="$SEUnwarpDir" \
  --t1samplespacing="$T1wSampleSpacing" \
  --t2samplespacing="$T2wSampleSpacing" \
  --unwarpdir="$UnwarpDir" \
  --gdcoeffs="$GradientDistortionCoeffs" \
  --avgrdcmethod="$AvgrdcSTRING" \
  --topupconfig="$TopupConfig" \
  --processing-mode="$ProcessingMode" \
  --printcom="$PRINTCOM" > "${QA_DIR}/PreFreeSurfer.txt"

SubjectDIR="${StudyFolder}/${Subject}/T1w"
T1wImage="${SubjectDIR}/T1w_acpc_dc_restore.nii.gz"
T1wImageBrain="${SubjectDIR}/T1w_acpc_dc_restore_brain.nii.gz"
if [ "$T2wInputImages" = "NONE" ]; then
  T2wImage="NONE"
else
  T2wImage="${SubjectDIR}/T2w_acpc_dc_restore.nii.gz"
fi

echo "Running FreeSurferPipeline for ${Subject}"
"${HCPPIPEDIR}/FreeSurfer/FreeSurferPipeline.sh" \
  --subject="$Subject" \
  --subjectDIR="$SubjectDIR" \
  --t1="$T1wImage" \
  --t1brain="$T1wImageBrain" \
  --t2="$T2wImage" \
  --processing-mode="$ProcessingMode" > "${QA_DIR}/FreeSurfer.txt"

SurfaceAtlasDIR="${HCPPIPEDIR_Templates}/standard_mesh_atlases"
GrayordinatesSpaceDIR="${HCPPIPEDIR_Templates}/91282_Greyordinates"
GrayordinatesResolutions="2"
HighResMesh="164"
LowResMeshes="32"
SubcorticalGrayLabels="${HCPPIPEDIR_Config}/FreeSurferSubcorticalLabelTableLut.txt"
FreeSurferLabels="${HCPPIPEDIR_Config}/FreeSurferAllLut.txt"
ReferenceMyelinMaps="${HCPPIPEDIR_Templates}/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii"
RegName="MSMSulc"

echo "Running PostFreeSurferPipeline for ${Subject}"
"${HCPPIPEDIR}/PostFreeSurfer/PostFreeSurferPipeline.sh" \
  --path="$StudyFolder" \
  --subject="$Subject" \
  --surfatlasdir="$SurfaceAtlasDIR" \
  --grayordinatesdir="$GrayordinatesSpaceDIR" \
  --grayordinatesres="$GrayordinatesResolutions" \
  --hiresmesh="$HighResMesh" \
  --lowresmesh="$LowResMeshes" \
  --subcortgraylabels="$SubcorticalGrayLabels" \
  --freesurferlabels="$FreeSurferLabels" \
  --refmyelinmaps="$ReferenceMyelinMaps" \
  --regname="$RegName" \
  --processing-mode="$ProcessingMode" > "${QA_DIR}/PostFreeSurfer.txt"

mkdir -p "$ANAT_BASE"
rm -rf "${ANAT_BASE}/T1w" "${ANAT_BASE}/T2w" "${ANAT_BASE}/MNINonLinear" "${ANAT_BASE}/qa"
[ -d "${StudyFolder}/${Subject}/T1w" ] && mv "${StudyFolder}/${Subject}/T1w" "$ANAT_BASE/"
[ -d "${StudyFolder}/${Subject}/T2w" ] && mv "${StudyFolder}/${Subject}/T2w" "$ANAT_BASE/"
[ -d "${StudyFolder}/${Subject}/MNINonLinear" ] && mv "${StudyFolder}/${Subject}/MNINonLinear" "$ANAT_BASE/"
[ -d "${QA_DIR}" ] && mv "${QA_DIR}" "${ANAT_BASE}/qa"

echo "HCP anatomical pipeline complete for ${Subject}"

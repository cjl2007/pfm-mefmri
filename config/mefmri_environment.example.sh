#!/usr/bin/env bash
# Copy into your study config folder and set EnvironmentScript in mefmri.sh to it.
# Load your site's environment modules first, or uncomment/edit the paths below.
# No installation is performed by this file.
# export FSLDIR="/path/to/fsl"
# export FREESURFER_HOME="/path/to/freesurfer"
# export FS_LICENSE="/path/to/license.txt"
# export PATH="/path/to/workbench/bin_linux64:$PATH"
# export MSMBINDIR="/path/to/directory/containing/msm"
# export PATH="/path/to/SimNIBS/bin:$PATH"  # provides charm if CHARM_BIN is empty

: "${MEDIR:?Load this environment through a configured pfm-mefmri study config}"
: "${FSLDIR:?Set FSLDIR here or load your site FSL environment}"
: "${FREESURFER_HOME:?Set FREESURFER_HOME here or load your site FreeSurfer environment}"
source "$FSLDIR/etc/fslconf/fsl.sh"
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
# The bundled HCP setup resolves Workbench from PATH and respects MSMBINDIR.
source "$MEDIR/HCPpipelines-master/Examples/Scripts/SetUpHCPPipeline.sh"

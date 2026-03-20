#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)

MEDIR=$1
Subject=$2
StudyFolder=$3

# Normalize MEDIR and StudyFolder to absolute paths.
if [[ "${MEDIR:0:1}" != "/" ]]; then
	MEDIR="$(cd "$MEDIR" && pwd)"
fi

# Normalize to absolute path so later `cd` calls cannot break relative paths.
if [[ "${StudyFolder:0:1}" != "/" ]]; then
	StudyFolder="$(cd "$StudyFolder" && pwd)"
fi

Subdir="$StudyFolder"/"$Subject"
SUBJECTS_DIR="$Subdir"/anat/T1w/ # note: this is used for "bbregister" calls;
AtlasTemplate=$4
if [[ "${AtlasTemplate:0:1}" != "/" ]]; then
	AtlasTemplate="$(cd "$(dirname "$AtlasTemplate")" && pwd)/$(basename "$AtlasTemplate")"
fi
DOF=$5
NTHREADS=$6
StartSession=$7
AtlasSpace=${8:-${AtlasSpace:-T1w}}
FuncDirName=${9:-${FUNC_DIRNAME:-rest}}
FuncFilePrefix=${10:-${FUNC_FILE_PREFIX:-Rest}}
ApplyN4Bias=${APPLY_N4_BIAS:-0}

case "${AtlasSpace}" in
	T1w|MNINonlinear) ;;
	*)
		echo "ERROR: mefmri_func_coreg.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)"
		exit 2
		;;
esac
echo "[coreg] AtlasSpace=${AtlasSpace} (cortical ribbon mask will be generated in selected atlas space)"
echo "[coreg] Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"

# Read the JSON sidecars for each scan and write the text files used later in preprocessing.

python3 "$MEDIR"/lib/find_epi_params.py \
--subdir "$Subdir" --func-name "$FuncDirName" --start-session "$StartSession"

# Create SBrefs (average of the first few echoes) for each scan.
# These serve as intermediate coregistration targets when needed.

# Create a working directory for SBref generation.
mkdir -p "$Subdir"/func/"$FuncDirName"/AverageSBref
WDIR="$Subdir"/func/"$FuncDirName"/AverageSBref

# count the number of sessions
sessions=("$Subdir"/func/unprocessed/"$FuncDirName"/session_*)
sessions=$(seq 1 1 "${#sessions[@]}")

# Iterate through sessions.
for s in $sessions ; do

	# count number of runs for this session;
	runs=("$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_*)
	runs=$(seq 1 1 "${#runs[@]}")

	# Iterate over runs.
	for r in $runs ; do 

		# Read echo times.
		te=$(cat "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/TE.txt)
		n_te=0

		# Iterate over echoes.
		for i in $te ; do

			# Track the current echo index.
			n_te=`expr $n_te + 1` 

			# If there is no single-band reference image, drop the initial non-steady-state volumes and create one.
			if [ ! -f "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref_S"$s"_R"$r"_E"$n_te".nii.gz ]; then
				fslroi "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/"$FuncFilePrefix"_S"$s"_R"$r"_E"$n_te".nii.gz "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref_S"$s"_R"$r"_E"$n_te".nii.gz 10 1
				echo 10 > "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/rmVols.txt
			fi

		done

		# Use the first echo to estimate the bias field.
		cp "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref*_E1.nii.gz "$WDIR"/TMP_1.nii.gz
		
		# Estimate field inhomogeneity and resample the bias field image (ANTs -> FSL orientation).
		if [[ "$ApplyN4Bias" -eq 1 ]]; then
			N4BiasFieldCorrection -d 3 -i "$WDIR"/TMP_1.nii.gz -o ["$WDIR"/TMP_restored.nii.gz, "$WDIR"/Bias_field_"$s"_"$r".nii.gz]
			flirt -in "$WDIR"/Bias_field_"$s"_"$r".nii.gz -ref "$WDIR"/TMP_1.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat -out "$WDIR"/Bias_field_"$s"_"$r".nii.gz -interp spline #
		fi

		# Reset the echo counter.
		n_te=0 

		# Iterate over echoes.
		for i in $te ; do

			# Skip longer echo times for the SBref average.
			if [[ $i < 60 ]] ; then 

				n_te=`expr $n_te + 1`
				cp "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref*_E"$n_te.nii".gz "$WDIR"/TMP_"$n_te".nii.gz
				if [[ "$ApplyN4Bias" -eq 1 ]]; then
					fslmaths "$WDIR"/TMP_"$n_te".nii.gz -div "$WDIR"/Bias_field_"$s"_"$r".nii.gz "$WDIR"/TMP_"$n_te".nii.gz # apply correction;
				fi

			fi

		done

		# Combine and average across echoes.
		fslmerge -t "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz "$WDIR"/TMP_*.nii.gz ##> /dev/null 2>&1
		fslmaths "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz -Tmean "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz
		cp "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz "$WDIR"/SBref_"$s"_"$r".nii.gz
		rm "$WDIR"/TMP* # remove helper files

	done

done

# Coregister all SBrefs and create an average SBref for cross-scan alignment.

# build a list of all SBrefs;
images=("$WDIR"/SBref_*.nii.gz)

# Count images and average if needed.
if [ "${#images[@]}" \> 1 ]; then

	# Align and average the single-band reference (SBref) images.
	"$MEDIR"/res0urces/FuncAverage -n -o "$Subdir"/func/xfms/rest/AvgSBref.nii.gz \
	"$WDIR"/SBref_*.nii.gz ##> /dev/null 2>&1 

else

	# copy over the lone single-band reference (SBref) image;
	cp "${images[0]}" "$Subdir"/func/xfms/rest/AvgSBref.nii.gz #> /dev/null 2>&1

fi

# Create a clean temporary copy of the FreeSurfer folder.
rm -rf "$Subdir"/anat/T1w/freesurfer #> /dev/null 2>&1
cp -rf "$Subdir"/anat/T1w/"$Subject" "$Subdir"/anat/T1w/freesurfer #> /dev/null 2>&1

# ------------------------------------------------------------------------------
# Ensure white.deformed exists (for bbregister --s freesurfer --surf white.deformed)
# This assumes your FS temp subject is flattened at: $Subdir/anat/T1w/freesurfer/{mri,surf,...}
# and that SUBJECTS_DIR is: $Subdir/anat/T1w  (parent containing the "freesurfer" subject folder)
# ------------------------------------------------------------------------------

export SUBJECTS_DIR="$Subdir/anat/T1w"
FS_SUBJ="freesurfer"

T1ACPC="$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz"
ORIGMGZ="$SUBJECTS_DIR/$FS_SUBJ/mri/orig.mgz"

LH_DEF="$SUBJECTS_DIR/$FS_SUBJ/surf/lh.white.deformed"
RH_DEF="$SUBJECTS_DIR/$FS_SUBJ/surf/rh.white.deformed"

if [ ! -f "$LH_DEF" ] || [ ! -f "$RH_DEF" ]; then
	echo "[INFO] white.deformed surfaces not found. Creating from FS white using header-based mapping to ACPC T1."

  # sanity checks
  if [ ! -f "$ORIGMGZ" ]; then
  	echo "[ERROR] Missing $ORIGMGZ"
  	exit 1
  fi
  if [ ! -f "$T1ACPC" ]; then
  	echo "[ERROR] Missing $T1ACPC"
  	exit 1
  fi
  if [ ! -f "$SUBJECTS_DIR/$FS_SUBJ/surf/lh.white" ] || [ ! -f "$SUBJECTS_DIR/$FS_SUBJ/surf/rh.white" ]; then
  	echo "[ERROR] Missing base white surfaces in $SUBJECTS_DIR/$FS_SUBJ/surf/"
  	exit 1
  fi

  REG_TMP="$(mktemp -p /tmp "${FS_SUBJ}_orig2acpc_XXXXXX.dat")"

  tkregister2 \
  --mov  "$T1ACPC" \
  --targ "$ORIGMGZ" \
  --noedit --regheader \
  --reg "$REG_TMP" || { echo "[ERROR] tkregister2 failed"; rm -f "$REG_TMP"; exit 1; }

   # Create lh.white.deformed / rh.white.deformed (coords expressed in T1ACPC geometry)
   mri_surf2surf \
   --s "$FS_SUBJ" --hemi lh \
   --sval-xyz white \
   --reg "$REG_TMP" \
   --tval-xyz "$T1ACPC" \
   --tval "$LH_DEF" || { echo "[ERROR] mri_surf2surf lh failed"; rm -f "$REG_TMP"; exit 1; }

   mri_surf2surf \
   --s "$FS_SUBJ" --hemi rh \
   --sval-xyz white \
   --reg "$REG_TMP" \
   --tval-xyz "$T1ACPC" \
   --tval "$RH_DEF" || { echo "[ERROR] mri_surf2surf rh failed"; rm -f "$REG_TMP"; exit 1; }

   rm -f "$REG_TMP"

  # verify outputs exist (this is what bbregister is looking for)
  if [ ! -f "$LH_DEF" ] || [ ! -f "$RH_DEF" ]; then
  	echo "[ERROR] Expected outputs not created:"
  	echo "  $LH_DEF"
  	echo "  $RH_DEF"
  	ls -la "$SUBJECTS_DIR/$FS_SUBJ/surf/" | head -n 50
  	exit 1
  fi

  echo "[INFO] Created white.deformed surfaces:"
  ls -la "$LH_DEF" "$RH_DEF"
else
	echo "[INFO] Found existing white.deformed surfaces; skipping."
fi

# define the effective echo spacing;
EchoSpacing=$(cat $Subdir/func/xfms/rest/EffectiveEchoSpacing.txt) 
AvgPEDIR=$(find "$Subdir"/func/"$FuncDirName" -type f -name PE.txt | sort | head -n 1 | xargs cat 2>/dev/null)
if [ -z "$AvgPEDIR" ]; then
	AvgPEDIR=${EPIREG_PEDIR:--y}
fi

# register average SBref image to T1-weighted anatomical image using FSL's EpiReg (correct for spatial distortions using average field map); 
"$MEDIR"/res0urces/epi_reg_dof --dof="$DOF" --epi="$Subdir"/func/xfms/rest/AvgSBref.nii.gz --t1="$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz --t1brain="$Subdir"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz --out="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg --fmap="$Subdir"/func/field_maps/Avg_FM_rads_acpc.nii.gz --fmapmag="$Subdir"/func/field_maps/Avg_FM_mag_acpc.nii.gz --fmapmagbrain="$Subdir"/func/field_maps/Avg_FM_mag_acpc_brain.nii.gz --echospacing="$EchoSpacing" --wmseg="$Subdir"/anat/T1w/"$Subject"/mri/white.nii.gz --nofmapreg --pedir="$AvgPEDIR" ##> /dev/null 2>&1
applywarp --interp=spline --in="$Subdir"/func/xfms/rest/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg.nii.gz --warp="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg_warp.nii.gz

# use BBRegister (BBR) to fine-tune the existing co-registration & output FSL style transformation matrix;
bbregister --s freesurfer --mov "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg.nii.gz --init-reg "$MEDIR"/res0urces/eye.dat --surf white.deformed --bold --reg "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.dat --6 --o "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.nii.gz ##> /dev/null 2>&1 
tkregister2 --s freesurfer --noedit --reg "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.dat --mov "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg.nii.gz --targ "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz --fslregout "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.mat ##> /dev/null 2>&1 

# add BBR step as post warp linear transformation & generate inverse warp;
convertwarp --warp1="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg_warp.nii.gz --postmat="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.mat --ref="$AtlasTemplate" --out="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz
applywarp --interp=spline --in="$Subdir"/func/xfms/rest/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR.nii.gz --warp="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz
invwarp --ref="$Subdir"/func/xfms/rest/AvgSBref.nii.gz -w "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz -o "$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_inv_warp.nii.gz # invert func --> T1w anatomical warp; includ. dc.;

# combine warps (distorted SBref image --> T1w_acpc & anatomical image in acpc --> MNI atlas)
convertwarp --ref="$AtlasTemplate" --warp1="$Subdir"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz --warp2="$Subdir"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz --out="$Subdir"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz
applywarp --interp=spline --in="$Subdir"/func/xfms/rest/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$Subdir"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR.nii.gz --warp="$Subdir"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz
invwarp -w "$Subdir"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz -o "$Subdir"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR_inv_warp.nii.gz --ref="$Subdir"/func/xfms/rest/AvgSBref.nii.gz # generate an inverse warp; atlas --> distorted SBref image 

# Also coregister individual SBrefs to the target anatomical image.
# This supports comparison of the average field map and scan-specific field maps.

# create & define the "CoregQA" folder;
mkdir -p "$Subdir"/func/qa/CoregQA #> /dev/null 2>&1

# count the number of sessions
Sessions=("$Subdir"/func/"$FuncDirName"/session_*)
Sessions=$(seq $StartSession 1 "${#sessions[@]}")

func () {

	# count number of runs for this session;
	runs=("$2"/func/"$FuncDirName"/session_"$6"/run_*)
	runs=$(seq 1 1 "${#runs[@]}")

	# Iterate over runs.
	for r in $runs ; do
		RunPEDIR=$(cat "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/PE.txt 2>/dev/null)
		if [ -z "$RunPEDIR" ]; then
			RunPEDIR="$7"
		fi

		# check to see if this scan has a field map or not;
		if [ -f "$2/func/field_maps/AllFMs/FM_rads_acpc_S"$6"_R"$r".nii.gz" ]; then

			# define the effective echo spacing;
			EchoSpacing=$(cat "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/EffectiveEchoSpacing.txt)
		
			# register average SBref image to T1-weighted anatomical image using FSL's EpiReg (correct for spatial distortions using scan-specific field map); 
			"$1"/res0urces/epi_reg_dof --dof="$4" --epi="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --t1="$2"/anat/T1w/T1w_acpc_dc_restore.nii.gz --t1brain="$2"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz --out="$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r" --fmap="$2"/func/field_maps/AllFMs/FM_rads_acpc_S"$6"_R"$r".nii.gz --fmapmag="$2"/func/field_maps/AllFMs/FM_mag_acpc_S"$6"_R"$r".nii.gz --fmapmagbrain="$2"/func/field_maps/AllFMs/FM_mag_acpc_brain_S"$6"_R"$r".nii.gz --echospacing="$EchoSpacing" --wmseg="$2"/anat/T1w/"$3"/mri/white.nii.gz --nofmapreg --pedir="$RunPEDIR" ##> /dev/null 2>&1
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r"_warp.nii.gz

			# Use BBRegister (BBR) to refine the existing coregistration and write an FSL-style transform.
			bbregister --s freesurfer --mov "$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --init-reg "$1"/res0urces/eye.dat --surf white.deformed --bold --reg "$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".dat --6 --o "$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz ##> /dev/null 2>&1 
			tkregister2 --s freesurfer --noedit --reg "$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".dat --mov "$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --targ "$2"/anat/T1w/T1w_acpc_dc_restore.nii.gz --fslregout "$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".mat ##> /dev/null 2>&1 

			# add BBR step as post warp linear transformation & generate inverse warp;
			convertwarp --warp1="$2"/func/xfms/rest/SBref2acpc_EpiReg_S"$6"_R"$r"_warp.nii.gz --postmat="$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".mat --ref="$5" --out="$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			mv "$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz "$2"/func/qa/CoregQA/SBref2acpc_EpiReg+BBR_ScanSpecificFM_S"$6"_R"$r".nii.gz
			
			# warp SBref image into MNI atlas volume space in a single spline warp; can be used for CoregQA
			convertwarp --ref="$5" --warp1="$2"/func/xfms/rest/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz --warp2="$2"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz --out="$2"/func/xfms/rest/SBref2nonlin_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$2"/func/qa/CoregQA/SBref2nonlin_EpiReg+BBR_ScanSpecificFM_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/rest/SBref2nonlin_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz

		fi

        # repeat warps (ACPC, MNI) but this time with the native --> acpc co-registration using an average field map;
        flirt -dof "$4" -in "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz -ref "$2"/func/xfms/rest/AvgSBref.nii.gz -out "$2"/func/qa/CoregQA/SBref2AvgSBref_S"$6"_R"$r".nii.gz -omat "$2"/func/qa/CoregQA/SBref2AvgSBref_S"$6"_R"$r".mat
        applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --premat="$2"/func/qa/CoregQA/SBref2AvgSBref_S"$6"_R"$r".mat --warp="$2"/func/xfms/rest/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz --out="$2"/func/qa/CoregQA/SBref2acpc_EpiReg+BBR_AvgFM_S"$6"_R"$r".nii.gz --ref="$5"
        applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --premat="$2"/func/qa/CoregQA/SBref2AvgSBref_S"$6"_R"$r".mat --warp="$2"/func/xfms/rest/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz --out="$2"/func/qa/CoregQA/SBref2nonlin_EpiReg+BBR_AvgFM_S"$6"_R"$r".nii.gz --ref="$5"

	done
}

export FuncDirName
export -f func # also coregister individual SBrefs to the target anatomical image
parallel --jobs $NTHREADS func ::: $MEDIR ::: $Subdir ::: $Subject ::: $DOF ::: $AtlasTemplate ::: $Sessions ::: "$AvgPEDIR" ##> /dev/null 2>&1  

# Write run-level pointer files needed by later stages.
# (brain mask and subcortical mask in functional space)

# T2w anatomicals (whole brain & brain extracted in Atlas Template space)
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T2w_acpc_dc_restore.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T2w_acpc_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T2w_acpc_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat

# T1w anatomicals (whole brain & brain extracted in Atlas Template space)
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T1w_acpc_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T1w_acpc_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$Subdir"/func/xfms/rest/T1w_acpc_brain_func.nii.gz -bin "$Subdir"/func/xfms/rest/T1w_acpc_brain_func_mask.nii.gz

# Resample high-quality cortical ribbon mask into functional atlas grid.
# Expected source: anat/T1w/CorticalRibbon.nii.gz (fallback to .ni.gz typo if present).
CorticalRibbonSrc="$Subdir"/anat/T1w/CorticalRibbon.nii.gz
if [ ! -f "$CorticalRibbonSrc" ] && [ -f "$Subdir"/anat/T1w/CorticalRibbon.ni.gz ]; then
	CorticalRibbonSrc="$Subdir"/anat/T1w/CorticalRibbon.ni.gz
fi
if [ ! -f "$CorticalRibbonSrc" ]; then
	AutoLhRibbon="$Subdir"/anat/T1w/"$Subject"/mri/lh.ribbon.mgz
	AutoRhRibbon="$Subdir"/anat/T1w/"$Subject"/mri/rh.ribbon.mgz
	if [ -f "$AutoLhRibbon" ] && [ -f "$AutoRhRibbon" ]; then
		echo "[WARN] CorticalRibbon.nii.gz missing; auto-building from FreeSurfer lh/rh.ribbon.mgz"
		TmpLhRibbon="$Subdir"/func/xfms/rest/.tmp_lh.ribbon.nii.gz
		TmpRhRibbon="$Subdir"/func/xfms/rest/.tmp_rh.ribbon.nii.gz
		CorticalRibbonSrc="$Subdir"/func/xfms/rest/.tmp_CorticalRibbon_auto.nii.gz
		mri_convert -i "$AutoLhRibbon" -o "$TmpLhRibbon" --like "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz
		mri_convert -i "$AutoRhRibbon" -o "$TmpRhRibbon" --like "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz
		fslmaths "$TmpLhRibbon" -add "$TmpRhRibbon" "$CorticalRibbonSrc"
		fslmaths "$CorticalRibbonSrc" -bin "$CorticalRibbonSrc"
		rm -f "$TmpLhRibbon" "$TmpRhRibbon"
	else
		echo "[WARN] CorticalRibbon source missing and FreeSurfer ribbon mgz unavailable; using T1w_acpc_brain_mask as fallback."
		CorticalRibbonSrc="$Subdir"/anat/T1w/T1w_acpc_brain_mask.nii.gz
	fi
fi
# T1w/ACPC functional-grid ribbon (legacy/default path).
flirt -interp nearestneighbour -in "$CorticalRibbonSrc" -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/CorticalRibbon_acpc_func_mask.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$Subdir"/func/xfms/rest/CorticalRibbon_acpc_func_mask.nii.gz -bin "$Subdir"/func/xfms/rest/CorticalRibbon_acpc_func_mask.nii.gz

# MNINonlinear functional-grid ribbon (uses nonlinear ACPC->standard warp).
NonlinWarp="$Subdir"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz
if [ -f "$NonlinWarp" ]; then
	applywarp --interp=nn --in="$CorticalRibbonSrc" --ref="$AtlasTemplate" --warp="$NonlinWarp" --out="$Subdir"/func/xfms/rest/CorticalRibbon_nonlin_func_mask.nii.gz
	fslmaths "$Subdir"/func/xfms/rest/CorticalRibbon_nonlin_func_mask.nii.gz -bin "$Subdir"/func/xfms/rest/CorticalRibbon_nonlin_func_mask.nii.gz
elif [[ "$AtlasSpace" == "MNINonlinear" ]]; then
	echo "ERROR: missing nonlinear warp required for AtlasSpace=MNINonlinear: $NonlinWarp"
	exit 2
else
	echo "[INFO] Nonlinear ribbon mask skipped (missing warp: $NonlinWarp)"
fi
rm -f "$Subdir"/func/xfms/rest/.tmp_CorticalRibbon_auto.nii.gz

# MNINonlinear anatomicals (whole brain & brain extracted in Atlas Template space)
if [ -f "$Subdir"/anat/MNINonLinear/T2w_restore.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T2w_restore.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T2w_nonlin_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
if [ -f "$Subdir"/anat/MNINonLinear/T2w_restore_brain.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T2w_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T2w_nonlin_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
if [ -f "$Subdir"/anat/MNINonLinear/T1w_restore.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T1w_restore.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T1w_nonlin_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T1w_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$Subdir"/func/xfms/rest/T1w_nonlin_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$Subdir"/func/xfms/rest/T1w_nonlin_brain_func.nii.gz -bin "$Subdir"/func/xfms/rest/T1w_nonlin_brain_func_mask.nii.gz

# Write run-level intermediate target and warp pointers in shell/Python-friendly text files.
# these are consumed by headmotion module.
ScanSpecificFM=${SCAN_SPECIFIC_FM:-}
if [[ -z "$ScanSpecificFM" ]]; then
	# Backward compatibility for older configs that still define COREG_POINTER_POLICY.
	# scan_specific_if_available -> 1, everything else -> 0.
	case "${COREG_POINTER_POLICY:-}" in
		scan_specific_if_available) ScanSpecificFM=1 ;;
		*) ScanSpecificFM=0 ;;
	esac
fi
if [[ "$ScanSpecificFM" != "0" && "$ScanSpecificFM" != "1" ]]; then
	echo "ERROR: SCAN_SPECIFIC_FM must be 0 or 1 (got '$ScanSpecificFM')"
	exit 2
fi
if [[ "$AtlasSpace" == "MNINonlinear" ]]; then
	PointerWarpSpace="nonlin"
else
	PointerWarpSpace="acpc"
fi

PointerLog="$Subdir/func/qa/CoregQA/CoregPointerSelection.tsv"
echo -e "session\trun\trho_avgfm\trho_scan_specific\tselection" > "$PointerLog"

for s in $Sessions ; do
	runs=("$Subdir"/func/"$FuncDirName"/session_"$s"/run_*)
	runs=$(seq 1 1 "${#runs[@]}")
	for r in $runs ; do
		AvgWarp="$Subdir/func/xfms/rest/AvgSBref2${PointerWarpSpace}_EpiReg+BBR_warp.nii.gz"
		ScanWarp="$Subdir/func/xfms/rest/SBref2${PointerWarpSpace}_EpiReg+BBR_S${s}_R${r}_warp.nii.gz"
		RunSBref="$Subdir/func/$FuncDirName/session_${s}/run_${r}/SBref.nii.gz"
		AvgSBref="$Subdir/func/xfms/rest/AvgSBref.nii.gz"
		TargetTxt="$Subdir/func/$FuncDirName/session_${s}/run_${r}/IntermediateCoregTarget.txt"
		WarpTxt="$Subdir/func/$FuncDirName/session_${s}/run_${r}/Intermediate2ACPCWarp.txt"

		if [[ "$ScanSpecificFM" == "1" && -f "$ScanWarp" ]]; then
			echo "$RunSBref" > "$TargetTxt"
			echo "$ScanWarp" > "$WarpTxt"
			echo -e "${s}\t${r}\tnan\tnan\tscan_specific" >> "$PointerLog"
		else
			echo "$AvgSBref" > "$TargetTxt"
			echo "$AvgWarp" > "$WarpTxt"
			echo -e "${s}\t${r}\tnan\tnan\tavgfm" >> "$PointerLog"
		fi
	done
done

# Remove the temporary FreeSurfer folder.
rm -rf "$Subdir"/anat/T1w/freesurfer/ 

# Generate subcortical ROIs in the ACPC/nonlinear functional grid.
"$MEDIR"/lib/make_precise_subcortical_labels.sh "$Subdir" "$AtlasTemplate" "$MEDIR"

echo "[INFO] Coreg module complete. MATLAB CoregQA post-steps are not used in this revised pipeline."

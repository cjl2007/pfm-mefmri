#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)

MEDIR=$1
Subject=$2
StudyFolder=$3
Subdir="$StudyFolder"/"$Subject"
AtlasTemplate=$4
DOF=$5
NTHREADS=$6
StartSession=$7
AtlasSpace=${8:-${AtlasSpace:-T1w}}
FuncDirName=${9:-${FUNC_DIRNAME:-rest}}
FuncFilePrefix=${10:-${FUNC_FILE_PREFIX:-Rest}}
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-rest}"

case "${AtlasSpace}" in
	T1w|MNINonlinear) ;;
	*)
		echo "ERROR: mefmri_func_headmotion.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)"
		exit 2
		;;
esac
echo "[headmotion] AtlasSpace=${AtlasSpace}"
echo "[headmotion] Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"

# count the number of sessions
sessions=("$Subdir"/func/"$FuncDirName"/session_*)
sessions=$(seq $StartSession 1 "${#sessions[@]}")

# sweep the sessions;
for s in $sessions ; do

	# count number of runs for this session;
	runs=("$Subdir"/func/"$FuncDirName"/session_"$s"/run_*)
	runs=$(seq 1 1 "${#runs[@]}")

	# Iterate over runs.
	for r in $runs ; do

		# "AllScans.txt" contains 
		# dir. paths to every scan. 
		echo /session_"$s"/run_"$r" \
		>> "$Subdir"/AllScans.txt  

	done

done

# define a list of directories;
AllScans=$(cat "$Subdir"/AllScans.txt) # note: this is used for parallel processing purposes.
rm "$Subdir"/AllScans.txt # remove helper file

func () {
	ApplyN4Bias=${APPLY_N4_BIAS:-0}
	KeepMCF=${HEADMOTION_KEEP_MCF:-0}
	local ProcDir="$3/func/$FuncDirName/$5"
	local UnprocDir="$3/func/unprocessed/$FuncDirName/$5"
	local BrainMask="$3/func/xfms/$FuncXfmsDir/T1w_acpc_brain_func_mask.nii.gz"
	if [[ "$AtlasSpace" == "MNINonlinear" ]]; then
		BrainMask="$3/func/xfms/$FuncXfmsDir/T1w_nonlin_brain_func_mask.nii.gz"
	elif [[ ! -f "$BrainMask" ]] && [[ -f "$3/func/xfms/$FuncXfmsDir/T1w_func_brain_mask.nii.gz" ]]; then
		# Compatibility fallback for older pre-refactor outputs.
		BrainMask="$3/func/xfms/$FuncXfmsDir/T1w_func_brain_mask.nii.gz"
	fi

	# Clean up motion-correction workspace from prior runs.
	rm -rf "$ProcDir"/MCF #> /dev/null 2>&1
	rm -rf "$ProcDir"/Rest_AVG_mcf.mat "$ProcDir"/Rest_AVG_mcf.mat+ # remove stale split/mat dirs from prior runs
	mkdir "$ProcDir"/vols/

	# Read acquisition parameters.
	te=$(cat "$ProcDir"/TE.txt)
	tr=$(cat "$ProcDir"/TR.txt)
	n_te=0

	# Read run-specific target and warp pointers.
	IntermediateCoregTarget=$(cat "$ProcDir"/IntermediateCoregTarget.txt)
	Intermediate2ACPCWarp=$(cat "$ProcDir"/Intermediate2ACPCWarp.txt)

	# Iterate over echoes.
	for i in $te ; do

		# Track the current echo index.
		n_te=`expr $n_te + 1`

		# skip the longer te;
		if [[ $i < 60 ]] ; then 

			# split original 4D resting-state file into single 3D vols.;
			fslsplit "$UnprocDir"/"$FuncFilePrefix"*_E"$n_te".nii.gz \
			"$ProcDir"/vols/E"$n_te"_

		fi
	
	done

	# Iterate over the individual volumes.
	for i in $(seq -f "%04g" 0 $((`fslnvols "$UnprocDir"/"$FuncFilePrefix"*_E1.nii.gz` - 1))) ; do

	  	# combine te;
	  	fslmerge -t "$ProcDir"/vols/AVG_"$i".nii.gz "$ProcDir"/vols/E*_"$i".nii.gz
		fslmaths "$ProcDir"/vols/AVG_"$i".nii.gz -Tmean "$ProcDir"/vols/AVG_"$i".nii.gz

	done

	# merge the images;
	fslmerge -t "$ProcDir"/Rest_AVG.nii.gz "$ProcDir"/vols/AVG_*.nii.gz # note: used for estimating head motion;
	fslmerge -t "$ProcDir"/Rest_E1.nii.gz "$ProcDir"/vols/E1_*.nii.gz # note: used for estimating (very rough)bias field;
	rm -rf "$ProcDir"/vols/ # remove temporary directory

	# Use the first echo to estimate the bias field.
	fslmaths "$ProcDir"/Rest_E1.nii.gz -Tmean "$ProcDir"/Mean.nii.gz
	if [[ "$ApplyN4Bias" -eq 1 ]]; then
		N4BiasFieldCorrection -d 3 -i "$ProcDir"/Mean.nii.gz -o ["$ProcDir"/Mean_Restored.nii.gz,"$ProcDir"/Bias_field.nii.gz] # estimate field inhomogeneity
	fi
	rm "$ProcDir"/Rest_E1.nii.gz # remove helper file

	if [[ "$ApplyN4Bias" -eq 1 ]]; then
		# resample bias field image (ANTs --> FSL orientation);
		flirt -in "$ProcDir"/Bias_field.nii.gz -ref "$ProcDir"/Mean.nii.gz -applyxfm \
		-init "$1"/res0urces/ident.mat -out "$ProcDir"/Bias_field.nii.gz -interp spline

		# remove signal bias;
		fslmaths "$ProcDir"/Rest_AVG.nii.gz \
		-div "$ProcDir"/Bias_field.nii.gz \
		"$ProcDir"/Rest_AVG.nii.gz
	fi

	# Remove helper files.
	rm -f "$ProcDir"/Mean*.nii.gz
	rm -f "$ProcDir"/Bias*.nii.gz

	# remove the first few volumes if needed;
	if [[ -f "$ProcDir"/rmVols.txt ]]; then
    	nVols=`fslnvols "$ProcDir"/Rest_AVG.nii.gz`
    	rmVols=$(cat "$ProcDir"/rmVols.txt)
		fslroi "$ProcDir"/Rest_AVG.nii.gz "$ProcDir"/Rest_AVG.nii.gz "$rmVols" `expr $nVols - $rmVols`
 	fi

	# Run an initial MCFLIRT pass to estimate motion before slice-time correction.
	mcflirt -dof "$4" -stages 3 -plots -in "$ProcDir"/Rest_AVG.nii.gz -r "$ProcDir"/SBref.nii.gz -out "$ProcDir"/MCF
	rm "$ProcDir"/MCF.nii.gz # remove .nii output; not used moving forward

	# perform slice time correction; using custom timing file;
	slicetimer -i "$ProcDir"/Rest_AVG.nii.gz \
	-o "$ProcDir"/Rest_AVG.nii.gz -r $tr \
	--tcustom="$ProcDir"/SliceTiming.txt

	# now run another MCFLIRT; specify average sbref as ref. vol & output transformation matrices;
	mcflirt -dof "$4" -mats -stages 3 -in "$ProcDir"/Rest_AVG.nii.gz -r "$IntermediateCoregTarget" -out "$ProcDir"/Rest_AVG_mcf
	rm "$ProcDir"/Rest_AVG*.nii.gz # delete intermediate images; not needed moving forward;

	# sweep all of the echoes; 
	for e in $(seq 1 1 "$n_te") ; do

		# copy over echo "e"; 
		cp "$UnprocDir"/"$FuncFilePrefix"*_E"$e".nii.gz \
		"$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz

		# remove the first few volumes if needed;
		if [[ -f "$ProcDir"/rmVols.txt ]]; then
			fslroi "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz \
			"$rmVols" `expr $nVols - $rmVols`
	 	fi

		# perform slice time correction using custom timing file;
		slicetimer -i "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz \
		--tcustom="$ProcDir"/SliceTiming.txt \
		-r $tr -o "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz

		# split original data into individual volumes;
		fslsplit "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz \
		"$ProcDir"/Rest_AVG_mcf.mat/vol_ -t

		# define affine transformation matrices and associated target images.
		mats=("$ProcDir"/Rest_AVG_mcf.mat/MAT_*)
	    images=("$ProcDir"/Rest_AVG_mcf.mat/vol_*.nii.gz)
		if [[ ${#images[@]} -ne ${#mats[@]} ]]; then
			echo "[ERROR] headmotion applywarp pairing mismatch for $5 echo $e: images=${#images[@]} mats=${#mats[@]}" 1>&2
			return 1
		fi

		# Warp image volumes into ACPC grid using legacy serial behavior.
		for (( i=0; i<${#images[@]}; i++ )); do
			applywarp --interp=spline --in="${images["$i"]}" --premat="${mats["$i"]}" \
			--warp="$Intermediate2ACPCWarp" --out="${images["$i"]}" --ref="$2"
		done

		# merge corrected images into a single file & perform a brain extraction
		fslmerge -t "$ProcDir"/"$FuncFilePrefix"_E"$e"_acpc.nii.gz "$ProcDir"/Rest_AVG_mcf.mat/*.nii.gz
		fslmaths "$ProcDir"/"$FuncFilePrefix"_E"$e"_acpc.nii.gz -mas "$BrainMask" \
		"$ProcDir"/"$FuncFilePrefix"_E"$e"_acpc.nii.gz # note: this step reduces file size, which is generally desirable but not absolutely needed.

		# Remove helper files.
		rm "$ProcDir"/Rest_AVG_mcf.mat/*.nii.gz # split volumes
		rm "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz # temporary copy of run data

	done

	# rename mcflirt transform dir.;
	rm -rf "$ProcDir"/MCF #> /dev/null 2>&1
	mv "$ProcDir"/*_mcf*.mat "$ProcDir"/MCF

	if [[ "$ApplyN4Bias" -eq 1 ]]; then
		# Use the first echo to estimate the bias field.
		fslmaths "$ProcDir"/"$FuncFilePrefix"_E1_acpc.nii.gz -Tmean "$ProcDir"/Mean.nii.gz
		fslmaths "$ProcDir"/Mean.nii.gz -thr 0 "$ProcDir"/Mean.nii.gz # remove any negative values introduced by spline interpolation;
		N4BiasFieldCorrection -d 3 -i "$ProcDir"/Mean.nii.gz -o ["$ProcDir"/Mean_Restored.nii.gz,"$ProcDir"/Bias_field.nii.gz] # estimate field inhomogeneity
		flirt -in "$ProcDir"/Bias_field.nii.gz -ref "$ProcDir"/Mean.nii.gz -applyxfm -init "$1"/res0urces/ident.mat -out "$ProcDir"/Bias_field.nii.gz -interp spline # resample bias field image (ANTs --> FSL orientation);

		# sweep all of the echoes;
		for e in $(seq 1 1 "$n_te") ; do

			# Correct signal inhomogeneity.
			fslmaths "$ProcDir"/"$FuncFilePrefix"_E"$e"_acpc.nii.gz \
			-div "$ProcDir"/Bias_field.nii.gz \
			"$ProcDir"/"$FuncFilePrefix"_E"$e"_acpc.nii.gz

		done
	fi

	# Remove helper files.
	rm -f "$ProcDir"/Mean*.nii.gz
	rm -f "$ProcDir"/Bias*.nii.gz

	# optional cleanup of motion transform workspace to avoid rerun clutter.
	if [[ "$KeepMCF" -eq 0 ]]; then
		rm -rf "$ProcDir"/MCF "$ProcDir"/Rest_AVG_mcf.mat "$ProcDir"/Rest_AVG_mcf.mat+
	fi

}

export FuncDirName FuncFilePrefix AtlasSpace
export -f func # correct for head motion and warp to atlas space in single spline warp
parallel --jobs $NTHREADS func ::: $MEDIR ::: $AtlasTemplate ::: $Subdir ::: $DOF ::: $AllScans #> /dev/null 2>&1  

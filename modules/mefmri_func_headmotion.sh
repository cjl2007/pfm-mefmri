#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)

set -euo pipefail
shopt -s nullglob

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
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-$FuncDirName}"
SBREF_ECHO_COMBINATION="${SBREF_ECHO_COMBINATION:-mean}"
SBREF_MAX_TE_MS="${SBREF_MAX_TE_MS:-60}"
SBREF_T2SMAP_MASK_THR="${SBREF_T2SMAP_MASK_THR:-1}"
SBREF_T2SMAP_THREADS="${SBREF_T2SMAP_THREADS:-1}"
SBREF_KEEP_INTERMEDIATES="${SBREF_KEEP_INTERMEDIATES:-0}"
TEDANA_ENV="${TEDANA_ENV:-mefmri_env}"
TEDANA_ACTIVATE_MODE="${TEDANA_ACTIVATE_MODE:-conda_activate}"
T2SMAP_BIN="${T2SMAP_BIN:-t2smap}"

case "${AtlasSpace}" in
	T1w|MNINonlinear) ;;
	*)
		echo "ERROR: mefmri_func_headmotion.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)"
		exit 2
		;;
esac
echo "[headmotion] AtlasSpace=${AtlasSpace}"
echo "[headmotion] Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"

float_lt() {
	python3 - "$1" "$2" <<'PY'
import sys
print("1" if float(sys.argv[1]) < float(sys.argv[2]) else "0")
PY
}

ms_to_s() {
	python3 - "$1" <<'PY'
import sys
print(f"{float(sys.argv[1]) / 1000.0:.8g}")
PY
}

resolve_tedana_bin() {
	local exe="$1"
	if [[ "$exe" == */* ]]; then
		[[ -x "$exe" ]] || { echo "ERROR: configured tedana executable is not executable: $exe" >&2; exit 2; }
		echo "$exe"
		return 0
	fi
	case "$TEDANA_ACTIVATE_MODE" in
		direct)
			command -v "$exe" || { echo "ERROR: could not find $exe on PATH with TEDANA_ACTIVATE_MODE=direct" >&2; exit 2; }
			;;
		conda_activate|conda_run)
			command -v conda >/dev/null 2>&1 || { echo "ERROR: conda is required for TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE" >&2; exit 2; }
			local base candidate
			base="$(conda info --base 2>/dev/null || true)"
			[[ -n "$base" ]] || { echo "ERROR: could not resolve conda base path" >&2; exit 2; }
			candidate="$base/envs/$TEDANA_ENV/bin/$exe"
			[[ -x "$candidate" ]] || { echo "ERROR: missing $exe in conda env '$TEDANA_ENV': $candidate" >&2; exit 2; }
			echo "$candidate"
			;;
		*)
			echo "ERROR: invalid TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE" >&2
			exit 2
			;;
	esac
}

mapfile -t sessions < <(
	find -L "$Subdir"/func/"$FuncDirName" -mindepth 1 -maxdepth 1 -type d -name 'session_*' -printf '%f\n' |
		sed 's/^session_//' |
		awk -v start="$StartSession" '$1 >= start' |
		sort -n
)

# sweep the sessions;
for s in "${sessions[@]}" ; do

	mapfile -t runs < <(
		find -L "$Subdir"/func/"$FuncDirName"/session_"$s" -mindepth 1 -maxdepth 1 -type d -name 'run_*' -printf '%f\n' |
			sed 's/^run_//' |
			sort -n
	)

	# Iterate over runs.
	for r in "${runs[@]}" ; do

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
	ReorientFuncToStd=${FUNC_REORIENT_TO_STD:-1}
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
	rm -rf "$ProcDir"/vols
	rm -f "$ProcDir"/Rest_AVG.nii.gz "$ProcDir"/Rest_AVG_mcf.nii.gz "$ProcDir"/Rest_E1.nii.gz
	rm -f "$ProcDir"/"$FuncFilePrefix"_E*.nii.gz "$ProcDir"/"$FuncFilePrefix"_E*_acpc.nii.gz
	mkdir "$ProcDir"/vols/

	reorient_to_std_inplace() {
		local img="$1"
		local tmp="${img%.nii.gz}_reorient_tmp.nii.gz"
		[[ "$ReorientFuncToStd" == "1" ]] || return 0
		fslreorient2std "$img" "$tmp"
		mv -f "$tmp" "$img"
	}

	combine_scout_echoes() {
		local method="$1"
		local out="$2"
		local tmp_dir="$3"
		shift 3
		local refs=("$@")
		[[ "${#refs[@]}" -gt 0 ]] || { echo "[ERROR] headmotion no echoes selected for scout combination" 1>&2; return 1; }
		case "$method" in
			mean)
				local sum="$tmp_dir/scout_sum.nii.gz"
				cp "${refs[0]}" "$sum"
				local idx
				for ((idx = 1; idx < ${#refs[@]}; idx++)); do
					fslmaths "$sum" -add "${refs[$idx]}" "$sum"
				done
				fslmaths "$sum" -div "${#refs[@]}" "$out"
				;;
			sos)
				local sum="$tmp_dir/scout_sos_sum.nii.gz"
				local sq idx
				fslmaths "${refs[0]}" -sqr "$sum"
				for ((idx = 1; idx < ${#refs[@]}; idx++)); do
					sq="$tmp_dir/scout_sos_sq_${idx}.nii.gz"
					fslmaths "${refs[$idx]}" -sqr "$sq"
					fslmaths "$sum" -add "$sq" "$sum"
				done
				fslmaths "$sum" -sqrt "$out"
				;;
			first)
				cp "${refs[0]}" "$out"
				;;
			last)
				cp "${refs[$((${#refs[@]} - 1))]}" "$out"
				;;
			t2s|paid)
				[[ "${#refs[@]}" -ge 2 ]] || { echo "[ERROR] headmotion $method scout combination requires at least two echoes" 1>&2; return 1; }
				[[ "${#selected_tes_seconds[@]}" -eq "${#refs[@]}" ]] || { echo "[ERROR] headmotion t2smap TE/reference count mismatch" 1>&2; return 1; }
				local t2smap_bin t2s_dir
				t2smap_bin="$(resolve_tedana_bin "$T2SMAP_BIN")"
				t2s_dir="$tmp_dir/t2smap_scout"
				mkdir -p "$t2s_dir" "$tmp_dir/mplconfig"
				fslmaths "${refs[0]}" -Tmean -thr "$SBREF_T2SMAP_MASK_THR" -bin "$tmp_dir/t2smap_scout_mask.nii.gz"
				MPLCONFIGDIR="$tmp_dir/mplconfig" "$t2smap_bin" \
					-d "${refs[@]}" \
					-e "${selected_tes_seconds[@]}" \
					--out-dir "$t2s_dir" \
					--prefix scout_ \
					--convention bids \
					--mask "$tmp_dir/t2smap_scout_mask.nii.gz" \
					--masktype none \
					--combmode "$method" \
					--n-threads "$SBREF_T2SMAP_THREADS" \
					--overwrite
				[[ -f "$t2s_dir/scout_desc-optcom_bold.nii.gz" ]] || { echo "[ERROR] t2smap did not create scout_desc-optcom_bold.nii.gz" 1>&2; return 1; }
				cp "$t2s_dir/scout_desc-optcom_bold.nii.gz" "$out"
				;;
			*)
				echo "[ERROR] unsupported SBREF_ECHO_COMBINATION=$method (use mean|sos|t2s|paid|first|last)" 1>&2
				return 1
				;;
		esac
	}

	# Read acquisition parameters.
	[[ -f "$ProcDir"/TE.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/TE.txt" 1>&2; return 1; }
	[[ -f "$ProcDir"/TR.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/TR.txt" 1>&2; return 1; }
	[[ -f "$ProcDir"/SliceTiming.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/SliceTiming.txt" 1>&2; return 1; }
	te=$(cat "$ProcDir"/TE.txt)
	tr=$(cat "$ProcDir"/TR.txt)
	n_te=0

	# Read run-specific target and warp pointers.
	[[ -f "$ProcDir"/IntermediateCoregTarget.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/IntermediateCoregTarget.txt" 1>&2; return 1; }
	IntermediateCoregTarget=$(cat "$ProcDir"/IntermediateCoregTarget.txt)
	if [[ "${FUNC_NOFIELDMAP_MODE:-0}" == "1" ]]; then
		[[ -f "$ProcDir"/Intermediate2ACPCMat.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/Intermediate2ACPCMat.txt" 1>&2; return 1; }
		Intermediate2ACPCMat=$(cat "$ProcDir"/Intermediate2ACPCMat.txt)
		[[ -f "$Intermediate2ACPCMat" ]] || { echo "[ERROR] headmotion missing affine target $Intermediate2ACPCMat" 1>&2; return 1; }
	else
		[[ -f "$ProcDir"/Intermediate2ACPCWarp.txt ]] || { echo "[ERROR] headmotion missing $ProcDir/Intermediate2ACPCWarp.txt" 1>&2; return 1; }
		Intermediate2ACPCWarp=$(cat "$ProcDir"/Intermediate2ACPCWarp.txt)
	fi

	# Copy selected echoes for the motion scout. In TOPUP/no-fieldmap mode this
	# scout intentionally remains in distorted/native EPI space.
	selected_refs=()
	selected_tes_seconds=()
	scout_sources="$ProcDir/HeadmotionScoutEchoSources.tsv"
	{
		echo -e "echo\tte_ms\tincluded\tsource_path\tworking_copy"
	} > "$scout_sources"
	for i in $te ; do

		# Track the current echo index.
		n_te=`expr $n_te + 1`

		raw_echo=( "$UnprocDir"/"$FuncFilePrefix"*_E"$n_te".nii.gz )
		[[ "${#raw_echo[@]}" -gt 0 ]] || { echo "[ERROR] headmotion missing raw echo $n_te in $UnprocDir" 1>&2; return 1; }
		include=1
		if [[ "$SBREF_MAX_TE_MS" != "none" ]]; then
			include="$(float_lt "$i" "$SBREF_MAX_TE_MS")"
		fi
		if [[ "$include" != "1" ]]; then
			echo -e "${n_te}\t${i}\t0\t${raw_echo[0]}\t" >> "$scout_sources"
			continue
		fi

		echo_copy="$ProcDir"/vols/E"$n_te"_source.nii.gz
		cp "${raw_echo[0]}" "$echo_copy"
		reorient_to_std_inplace "$echo_copy"
		selected_refs+=("$echo_copy")
		selected_tes_seconds+=("$(ms_to_s "$i")")
		echo -e "${n_te}\t${i}\t1\t${raw_echo[0]}\t${echo_copy}" >> "$scout_sources"
	done

	raw_e1=( "$UnprocDir"/"$FuncFilePrefix"*_E1.nii.gz )
	[[ "${#raw_e1[@]}" -gt 0 ]] || { echo "[ERROR] headmotion missing first-echo input in $UnprocDir" 1>&2; return 1; }
	combine_scout_echoes "$SBREF_ECHO_COMBINATION" "$ProcDir"/Rest_AVG.nii.gz "$ProcDir"/vols "${selected_refs[@]}"
	cp "$ProcDir"/vols/E1_source.nii.gz "$ProcDir"/Rest_E1.nii.gz

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
	if [[ "$SBREF_KEEP_INTERMEDIATES" == "1" ]]; then
		echo "$ProcDir/vols" > "$ProcDir/HeadmotionScoutIntermediatesDir.txt"
	else
		rm -rf "$ProcDir"/vols/ # remove temporary directory
	fi

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
		raw_echo=( "$UnprocDir"/"$FuncFilePrefix"*_E"$e".nii.gz )
		[[ "${#raw_echo[@]}" -gt 0 ]] || { echo "[ERROR] headmotion missing raw echo $e in $UnprocDir" 1>&2; return 1; }
		cp "${raw_echo[0]}" \
		"$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz
		reorient_to_std_inplace "$ProcDir"/"$FuncFilePrefix"_E"$e".nii.gz

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

		# Apply motion + atlas transforms volume-wise.
		for (( i=0; i<${#images[@]}; i++ )); do
			if [[ "${FUNC_NOFIELDMAP_MODE:-0}" == "1" ]]; then
				combined_mat="${mats["$i"]}_to_acpc.mat"
				convert_xfm -omat "$combined_mat" -concat "$Intermediate2ACPCMat" "${mats["$i"]}"
				flirt -in "${images["$i"]}" -ref "$2" -out "${images["$i"]}" \
				-applyxfm -init "$combined_mat" -interp spline
				rm -f "$combined_mat"
			else
				applywarp --interp=spline --in="${images["$i"]}" --premat="${mats["$i"]}" \
				--warp="$Intermediate2ACPCWarp" --out="${images["$i"]}" --ref="$2"
			fi
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

export FuncDirName FuncFilePrefix FuncXfmsDir AtlasSpace FUNC_NOFIELDMAP_MODE
export SBREF_ECHO_COMBINATION SBREF_MAX_TE_MS SBREF_T2SMAP_MASK_THR SBREF_T2SMAP_THREADS SBREF_KEEP_INTERMEDIATES
export TEDANA_ENV TEDANA_ACTIVATE_MODE T2SMAP_BIN
export -f func float_lt ms_to_s resolve_tedana_bin # correct for head motion and warp to atlas space in single spline warp
parallel --jobs $NTHREADS func ::: $MEDIR ::: $AtlasTemplate ::: $Subdir ::: $DOF ::: $AllScans #> /dev/null 2>&1  

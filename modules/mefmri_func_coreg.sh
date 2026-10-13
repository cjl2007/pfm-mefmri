#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)

set -euo pipefail
shopt -s nullglob

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
FuncXfmsDir=${FUNC_XFMS_DIRNAME:-${FuncDirName}}
ApplyN4Bias=${APPLY_N4_BIAS:-0}
COREG_PYTHON="${COREG_PYTHON:-${PIPELINE_PYTHON:-python3}}"
FUNC_NOFIELDMAP_MODE="${FUNC_NOFIELDMAP_MODE:-0}"
FM_OUT_DIR_REL="${FM_OUT_DIR_REL:-func/${FuncDirName}/field_maps}"
FM_DIR="$Subdir/$FM_OUT_DIR_REL"
QC_OUTPUT_DIRNAME="${QC_OUTPUT_DIRNAME:-QC}"
COREG_QA_DIR="${COREG_QA_DIR:-$Subdir/func/${FuncDirName}/${QC_OUTPUT_DIRNAME}/Coreg}"

case "${AtlasSpace}" in
	T1w|MNINonlinear) ;;
	*)
		echo "ERROR: mefmri_func_coreg.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)"
		exit 2
		;;
esac
echo "[coreg] AtlasSpace=${AtlasSpace} (cortical ribbon mask will be generated in selected atlas space)"
echo "[coreg] Functional naming: func/${FuncDirName}, prefix ${FuncFilePrefix}_*"
echo "[coreg] Transform naming: func/xfms/${FuncXfmsDir}"
echo "[coreg] Fieldmap directory: ${FM_DIR}"
echo "[coreg] Coreg QA directory: ${COREG_QA_DIR}"
echo "[coreg] FUNC_NOFIELDMAP_MODE=${FUNC_NOFIELDMAP_MODE}"

ensure_white_wmseg() {
	local fs_mri_dir="$Subdir/anat/T1w/$Subject/mri"
	local aseg_mgz="$fs_mri_dir/aparc+aseg.mgz"
	local white_mgz="$fs_mri_dir/white.mgz"
	local white_nii="$fs_mri_dir/white.nii.gz"
	if [[ -f "$white_nii" ]]; then
		echo "$white_nii"
		return 0
	fi
	[[ -f "$aseg_mgz" ]] || {
		echo "ERROR: missing FreeSurfer aseg for WM segmentation: $aseg_mgz" >&2
		exit 1
	}
	mri_binarize --i "$aseg_mgz" --wm --o "$white_mgz" >/dev/null 2>&1
	mri_convert -i "$white_mgz" -o "$white_nii" --like "$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz" >/dev/null 2>&1
	[[ -f "$white_nii" ]] || {
		echo "ERROR: failed to create WM segmentation: $white_nii" >&2
		exit 1
	}
	echo "$white_nii"
}

# Read the JSON sidecars for each scan and write the text files used later in preprocessing.

"$COREG_PYTHON" "$MEDIR"/lib/find_epi_params.py \
--subdir "$Subdir" --func-name "$FuncDirName" --func-prefix "$FuncFilePrefix" --start-session "$StartSession" --qc-dirname "$QC_OUTPUT_DIRNAME" \
$( [[ "$FUNC_NOFIELDMAP_MODE" == "1" ]] && printf '%s' "--no-fieldmap-mode" )
WMSEG_NII="$(ensure_white_wmseg)"
SBREF_FALLBACK_SKIP_TRS="${SBREF_FALLBACK_SKIP_TRS:-10}"
SBREF_REORIENT_TO_STD="${SBREF_REORIENT_TO_STD:-1}"
SBREF_ECHO_COMBINATION="${SBREF_ECHO_COMBINATION:-${MEDIC_SBREF_COMBINATION:-mean}}"
SBREF_MAX_TE_MS="${SBREF_MAX_TE_MS:-${MEDIC_SBREF_MAX_TE_MS:-60}}"
SBREF_T2SMAP_MASK_THR="${SBREF_T2SMAP_MASK_THR:-${MEDIC_SBREF_T2SMAP_MASK_THR:-1}}"
SBREF_T2SMAP_THREADS="${SBREF_T2SMAP_THREADS:-${MEDIC_SBREF_T2SMAP_THREADS:-1}}"
SBREF_KEEP_INTERMEDIATES="${SBREF_KEEP_INTERMEDIATES:-${MEDIC_SBREF_KEEP_INTERMEDIATES:-0}}"
TEDANA_ENV="${TEDANA_ENV:-mefmri_env}"
TEDANA_ACTIVATE_MODE="${TEDANA_ACTIVATE_MODE:-conda_activate}"
T2SMAP_BIN="${T2SMAP_BIN:-t2smap}"

reorient_to_std_inplace() {
	local img="$1"
	local tmp="${img%.nii.gz}_reorient_tmp.nii.gz"
	[[ "$SBREF_REORIENT_TO_STD" == "1" ]] || return 0
	[[ -f "$img" ]] || return 1
	fslreorient2std "$img" "$tmp"
	mv -f "$tmp" "$img"
}

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
		[[ -x "$exe" ]] || {
			echo "ERROR: configured tedana executable is not executable: $exe" >&2
			exit 2
		}
		echo "$exe"
		return 0
	fi
	case "$TEDANA_ACTIVATE_MODE" in
		direct)
			command -v "$exe" || {
				echo "ERROR: could not find $exe on PATH with TEDANA_ACTIVATE_MODE=direct" >&2
				exit 2
			}
			;;
		conda_activate|conda_run)
			command -v conda >/dev/null 2>&1 || {
				echo "ERROR: conda is required for TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE" >&2
				exit 2
			}
			local base candidate
			base="$(conda info --base 2>/dev/null || true)"
			[[ -n "$base" ]] || {
				echo "ERROR: could not resolve conda base path" >&2
				exit 2
			}
			candidate="$base/envs/$TEDANA_ENV/bin/$exe"
			[[ -x "$candidate" ]] || {
				echo "ERROR: missing $exe in conda env '$TEDANA_ENV': $candidate" >&2
				exit 2
			}
			echo "$candidate"
			;;
		*)
			echo "ERROR: invalid TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE" >&2
			exit 2
			;;
	esac
}

make_4d_copy() {
	local py_bin="$1"
	local input="$2"
	local output="$3"
	"$py_bin" - "$input" "$output" <<'PY'
import sys
import nibabel as nb

img = nb.load(sys.argv[1])
data = img.get_fdata(dtype="float32")
if data.ndim == 3:
    data = data[..., None]
nb.save(nb.Nifti1Image(data, img.affine, img.header), sys.argv[2])
PY
}

combine_sbref_t2smap() {
	local combmode="$1"
	local out="$2"
	local tmp_dir="$3"
	shift 3
	local refs=("$@")

	[[ "${#COMBINE_TES_SECONDS[@]}" -eq "${#refs[@]}" ]] || {
		echo "ERROR: t2smap TE/reference count mismatch" >&2
		exit 2
	}
	[[ "${#refs[@]}" -ge 2 ]] || {
		echo "ERROR: t2smap SBRef combination requires at least two echoes" >&2
		exit 2
	}

	local t2smap_bin py_bin t2s_dir inputs=()
	t2smap_bin="$(resolve_tedana_bin "$T2SMAP_BIN")"
	py_bin="$(resolve_tedana_bin python)"
	t2s_dir="$tmp_dir/t2smap"
	mkdir -p "$t2s_dir" "$tmp_dir/mplconfig"

	fslmaths "${refs[0]}" -thr "$SBREF_T2SMAP_MASK_THR" -bin "$tmp_dir/t2smap_mask.nii.gz"

	local idx ref4
	for ((idx = 0; idx < ${#refs[@]}; idx++)); do
		ref4="$tmp_dir/t2smap_echo_$((idx + 1))_4d.nii.gz"
		make_4d_copy "$py_bin" "${refs[$idx]}" "$ref4"
		inputs+=("$ref4")
	done

	MPLCONFIGDIR="$tmp_dir/mplconfig" "$t2smap_bin" \
		-d "${inputs[@]}" \
		-e "${COMBINE_TES_SECONDS[@]}" \
		--out-dir "$t2s_dir" \
		--prefix sbref_ \
		--convention bids \
		--mask "$tmp_dir/t2smap_mask.nii.gz" \
		--masktype none \
		--combmode "$combmode" \
		--n-threads "$SBREF_T2SMAP_THREADS" \
		--overwrite

	[[ -f "$t2s_dir/sbref_desc-optcom_bold.nii.gz" ]] || {
		echo "ERROR: t2smap did not create sbref_desc-optcom_bold.nii.gz" >&2
		exit 2
	}
	fslroi "$t2s_dir/sbref_desc-optcom_bold.nii.gz" "$out" 0 1
}

combine_sbref_echoes() {
	local method="$1"
	local out="$2"
	local tmp_dir="$3"
	shift 3
	local refs=("$@")

	[[ "${#refs[@]}" -gt 0 ]] || {
		echo "ERROR: no SBRef echoes were selected for combination" >&2
		exit 2
	}

	case "$method" in
		mean)
			fslmerge -t "$tmp_dir/all_echo_refs.nii.gz" "${refs[@]}"
			fslmaths "$tmp_dir/all_echo_refs.nii.gz" -Tmean "$out"
			;;
		sos)
			local sum="$tmp_dir/sos_sum.nii.gz"
			local sq
			fslmaths "${refs[0]}" -sqr "$sum"
			local idx
			for ((idx = 1; idx < ${#refs[@]}; idx++)); do
				sq="$tmp_dir/sos_sq_${idx}.nii.gz"
				fslmaths "${refs[$idx]}" -sqr "$sq"
				fslmaths "$sum" -add "$sq" "$sum"
			done
			fslmaths "$sum" -sqrt "$out"
			;;
		t2s|paid)
			combine_sbref_t2smap "$method" "$out" "$tmp_dir" "${refs[@]}"
			;;
		first)
			cp "${refs[0]}" "$out"
			;;
		last)
			cp "${refs[$((${#refs[@]} - 1))]}" "$out"
			;;
		*)
			echo "ERROR: unsupported SBREF_ECHO_COMBINATION=$method (use mean|sos|t2s|paid|first|last)" >&2
			exit 2
			;;
	esac
}

# Create SBrefs (average of the first few echoes) for each scan.
# These serve as intermediate coregistration targets when needed.

# Create a working directory for SBref generation.
mkdir -p "$Subdir"/func/"$FuncDirName"/AverageSBref
mkdir -p "$Subdir"/func/xfms/"$FuncXfmsDir"
WDIR="$Subdir"/func/"$FuncDirName"/AverageSBref
XfmsDir="$Subdir"/func/xfms/"$FuncXfmsDir"

mapfile -t sessions < <(
	find -L "$Subdir"/func/unprocessed/"$FuncDirName" -mindepth 1 -maxdepth 1 -type d -name 'session_*' -printf '%f\n' |
		sed 's/^session_//' |
		sort -n
)

# Iterate through sessions.
for s in "${sessions[@]}" ; do

	mapfile -t runs < <(
		find -L "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s" -mindepth 1 -maxdepth 1 -type d -name 'run_*' -printf '%f\n' |
			sed 's/^run_//' |
			sort -n
	)

	# Iterate over runs.
	for r in "${runs[@]}" ; do

		# Read echo times.
		[[ -f "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/TE.txt ]] || {
			echo "ERROR: missing TE.txt for session_${s}/run_${r}. Run input validation/import before coreg." >&2
			exit 1
		}
		te=$(cat "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/TE.txt)
		n_te=0

		# Iterate over echoes.
		for i in $te ; do

			# Track the current echo index.
			n_te=`expr $n_te + 1` 

			# If there is no single-band reference image, drop the initial non-steady-state volumes and create one.
			if [ ! -f "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref_S"$s"_R"$r"_E"$n_te".nii.gz ]; then
				[[ -f "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/"$FuncFilePrefix"_S"$s"_R"$r"_E"$n_te".nii.gz ]] || {
					echo "ERROR: missing raw echo input for session_${s}/run_${r}/echo_${n_te} under func/unprocessed/${FuncDirName}. Import BIDS/raw inputs first." >&2
					exit 1
				}
				fslroi "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/"$FuncFilePrefix"_S"$s"_R"$r"_E"$n_te".nii.gz "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref_S"$s"_R"$r"_E"$n_te".nii.gz "$SBREF_FALLBACK_SKIP_TRS" 1
				echo "$SBREF_FALLBACK_SKIP_TRS" > "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/rmVols.txt
			fi
			reorient_to_std_inplace "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref_S"$s"_R"$r"_E"$n_te".nii.gz

		done

		# Use the first echo to estimate the bias field.
		sbref_e1=( "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref*_E1.nii.gz )
		[[ "${#sbref_e1[@]}" -gt 0 ]] || {
			echo "ERROR: failed to create/find SBref E1 for session_${s}/run_${r}." >&2
			exit 1
		}
		cp "${sbref_e1[0]}" "$WDIR"/TMP_1.nii.gz
		
		# Estimate field inhomogeneity and resample the bias field image (ANTs -> FSL orientation).
		if [[ "$ApplyN4Bias" -eq 1 ]]; then
			N4BiasFieldCorrection -d 3 -i "$WDIR"/TMP_1.nii.gz -o ["$WDIR"/TMP_restored.nii.gz, "$WDIR"/Bias_field_"$s"_"$r".nii.gz]
			flirt -in "$WDIR"/Bias_field_"$s"_"$r".nii.gz -ref "$WDIR"/TMP_1.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat -out "$WDIR"/Bias_field_"$s"_"$r".nii.gz -interp spline #
		fi

		# Combine selected echoes into the run-level SBref.
		n_te=0
		sbref_tmp="$WDIR"/sbref_"$s"_"$r"_tmp
		rm -rf "$sbref_tmp"
		mkdir -p "$sbref_tmp"
		selected_refs=()
		selected_tes_seconds=()
		sources_tsv="$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref_echo_sources.tsv
		{
			echo -e "echo\tte_ms\tincluded\tsource_type\tsource_path\tworking_copy"
		} > "$sources_tsv"

		for i in $te ; do
			n_te=`expr $n_te + 1`
			include=1
			if [[ "$SBREF_MAX_TE_MS" != "none" ]]; then
				include="$(float_lt "$i" "$SBREF_MAX_TE_MS")"
			fi

			if [[ "$include" != "1" ]]; then
				echo -e "${n_te}\t${i}\t0\tskipped_te_ge_${SBREF_MAX_TE_MS}\t\t" >> "$sources_tsv"
				continue
			fi

			sbref_echo=( "$Subdir"/func/unprocessed/"$FuncDirName"/session_"$s"/run_"$r"/SBref*_E"$n_te".nii.gz )
			[[ "${#sbref_echo[@]}" -gt 0 ]] || {
				echo "ERROR: missing SBref echo ${n_te} for session_${s}/run_${r}." >&2
				exit 1
			}
			one="$sbref_tmp"/echo_"$n_te".nii.gz
			cp "${sbref_echo[0]}" "$one"
			if [[ "$ApplyN4Bias" -eq 1 ]]; then
				fslmaths "$one" -div "$WDIR"/Bias_field_"$s"_"$r".nii.gz "$one"
			fi
			selected_refs+=("$one")
			selected_tes_seconds+=("$(ms_to_s "$i")")
			echo -e "${n_te}\t${i}\t1\tSBref\t${sbref_echo[0]}\t${one}" >> "$sources_tsv"
		done

		COMBINE_TES_SECONDS=("${selected_tes_seconds[@]}")
		combine_sbref_echoes "$SBREF_ECHO_COMBINATION" "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz "$sbref_tmp" "${selected_refs[@]}"
		reorient_to_std_inplace "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz
		cp "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref.nii.gz "$WDIR"/SBref_"$s"_"$r".nii.gz
		if [[ "$SBREF_KEEP_INTERMEDIATES" == "1" ]]; then
			echo "$sbref_tmp" > "$Subdir"/func/"$FuncDirName"/session_"$s"/run_"$r"/SBref_intermediates_dir.txt
		else
			rm -rf "$sbref_tmp"
		fi
		rm -f "$WDIR"/TMP* # remove helper files

	done

done

# Coregister all SBrefs and create an average SBref for cross-scan alignment.

# build a list of all SBrefs;
images=("$WDIR"/SBref_*.nii.gz)

# Count images and average if needed.
if [ "${#images[@]}" \> 1 ]; then

	# Align and average the single-band reference (SBref) images.
	"$MEDIR"/res0urces/FuncAverage -n -o "$XfmsDir"/AvgSBref.nii.gz \
	"$WDIR"/SBref_*.nii.gz ##> /dev/null 2>&1 

else

	# copy over the lone single-band reference (SBref) image;
	cp "${images[0]}" "$XfmsDir"/AvgSBref.nii.gz #> /dev/null 2>&1

fi
reorient_to_std_inplace "$XfmsDir"/AvgSBref.nii.gz

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

AvgPEDIR=$(find -L "$Subdir"/func/"$FuncDirName" -type f -name PE.txt | sort | head -n 1 | xargs cat 2>/dev/null)
if [ -z "$AvgPEDIR" ]; then
	AvgPEDIR=${EPIREG_PEDIR:--y}
fi

# register average SBref image to T1-weighted anatomical image.
if [[ "$FUNC_NOFIELDMAP_MODE" == "1" ]]; then
	"$MEDIR"/res0urces/epi_reg_dof --dof="$DOF" --epi="$XfmsDir"/AvgSBref.nii.gz --t1="$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz --t1brain="$Subdir"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz --out="$XfmsDir"/AvgSBref2acpc_EpiReg --wmseg="$WMSEG_NII" ##> /dev/null 2>&1
	convertwarp --ref="$AtlasTemplate" --premat="$XfmsDir"/AvgSBref2acpc_EpiReg.mat --out="$XfmsDir"/AvgSBref2acpc_EpiReg_warp.nii.gz
else
	[[ -f "$XfmsDir"/EffectiveEchoSpacing.txt ]] || {
		echo "ERROR: missing func/xfms/${FuncXfmsDir}/EffectiveEchoSpacing.txt. Re-run validation/import so JSON metadata are available." >&2
		exit 1
	}
	required_fmaps=(
		"$FM_DIR/Avg_FM_rads_acpc.nii.gz"
		"$FM_DIR/Avg_FM_mag_acpc.nii.gz"
		"$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz"
	)
	for fmap_file in "${required_fmaps[@]}"; do
		[[ -f "$fmap_file" ]] || {
			echo "ERROR: missing processed fieldmap: $fmap_file" >&2
			echo "Re-run the fieldmaps module before coreg." >&2
			exit 1
		}
	done
	EchoSpacing=$(cat "$XfmsDir"/EffectiveEchoSpacing.txt)
	"$MEDIR"/res0urces/epi_reg_dof --dof="$DOF" --epi="$XfmsDir"/AvgSBref.nii.gz --t1="$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz --t1brain="$Subdir"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz --out="$XfmsDir"/AvgSBref2acpc_EpiReg --fmap="$FM_DIR"/Avg_FM_rads_acpc.nii.gz --fmapmag="$FM_DIR"/Avg_FM_mag_acpc.nii.gz --fmapmagbrain="$FM_DIR"/Avg_FM_mag_acpc_brain.nii.gz --echospacing="$EchoSpacing" --wmseg="$WMSEG_NII" --nofmapreg --pedir="$AvgPEDIR" ##> /dev/null 2>&1
fi
applywarp --interp=spline --in="$XfmsDir"/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$XfmsDir"/AvgSBref2acpc_EpiReg.nii.gz --warp="$XfmsDir"/AvgSBref2acpc_EpiReg_warp.nii.gz

# use BBRegister (BBR) to fine-tune the existing co-registration & output FSL style transformation matrix;
bbregister --s freesurfer --mov "$XfmsDir"/AvgSBref2acpc_EpiReg.nii.gz --init-reg "$MEDIR"/res0urces/eye.dat --surf white.deformed --bold --reg "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.dat --6 --o "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.nii.gz ##> /dev/null 2>&1 
tkregister2 --s freesurfer --noedit --reg "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.dat --mov "$XfmsDir"/AvgSBref2acpc_EpiReg.nii.gz --targ "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz --fslregout "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.mat ##> /dev/null 2>&1 

# add BBR step as post warp linear transformation & generate inverse warp;
convertwarp --warp1="$XfmsDir"/AvgSBref2acpc_EpiReg_warp.nii.gz --postmat="$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.mat --ref="$AtlasTemplate" --out="$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz
applywarp --interp=spline --in="$XfmsDir"/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.nii.gz --warp="$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz
invwarp --ref="$XfmsDir"/AvgSBref.nii.gz -w "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz -o "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_inv_warp.nii.gz # invert func --> T1w anatomical warp; includ. dc.;
convert_xfm -omat "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_full.mat \
-concat "$XfmsDir"/AvgSBref2acpc_EpiReg+BBR.mat \
"$XfmsDir"/AvgSBref2acpc_EpiReg.mat

# combine warps (distorted SBref image --> T1w_acpc & anatomical image in acpc --> MNI atlas)
convertwarp --ref="$AtlasTemplate" --warp1="$XfmsDir"/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz --warp2="$Subdir"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz --out="$XfmsDir"/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz
applywarp --interp=spline --in="$XfmsDir"/AvgSBref.nii.gz --ref="$AtlasTemplate" --out="$XfmsDir"/AvgSBref2nonlin_EpiReg+BBR.nii.gz --warp="$XfmsDir"/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz
invwarp -w "$XfmsDir"/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz -o "$XfmsDir"/AvgSBref2nonlin_EpiReg+BBR_inv_warp.nii.gz --ref="$XfmsDir"/AvgSBref.nii.gz # generate an inverse warp; atlas --> distorted SBref image 

# Also coregister individual SBrefs to the target anatomical image.
# This supports comparison of the average field map and scan-specific field maps.

# create & define the task-specific coregistration QA folder;
mkdir -p "$COREG_QA_DIR" #> /dev/null 2>&1

mapfile -t Sessions < <(
	find -L "$Subdir"/func/"$FuncDirName" -mindepth 1 -maxdepth 1 -type d -name 'session_*' -printf '%f\n' |
		sed 's/^session_//' |
		awk -v start="$StartSession" '$1 >= start' |
		sort -n
)

func () {

	mapfile -t runs < <(
		find -L "$2"/func/"$FuncDirName"/session_"$6" -mindepth 1 -maxdepth 1 -type d -name 'run_*' -printf '%f\n' |
			sed 's/^run_//' |
			sort -n
	)

	# Iterate over runs.
	for r in "${runs[@]}" ; do
		RunPEDIR=$(cat "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/PE.txt 2>/dev/null)
		if [ -z "$RunPEDIR" ]; then
			RunPEDIR="$7"
		fi

		# check to see if this scan has a field map or not;
		if [ -f "$FM_DIR/AllFMs/FM_rads_acpc_S"$6"_R"$r".nii.gz" ]; then

			# define the effective echo spacing;
			EchoSpacing=$(cat "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/EffectiveEchoSpacing.txt)
		
			# register average SBref image to T1-weighted anatomical image using FSL's EpiReg (correct for spatial distortions using scan-specific field map); 
			"$1"/res0urces/epi_reg_dof --dof="$4" --epi="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --t1="$2"/anat/T1w/T1w_acpc_dc_restore.nii.gz --t1brain="$2"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz --out="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r" --fmap="$FM_DIR"/AllFMs/FM_rads_acpc_S"$6"_R"$r".nii.gz --fmapmag="$FM_DIR"/AllFMs/FM_mag_acpc_S"$6"_R"$r".nii.gz --fmapmagbrain="$FM_DIR"/AllFMs/FM_mag_acpc_brain_S"$6"_R"$r".nii.gz --echospacing="$EchoSpacing" --wmseg="$WMSEG_NII" --nofmapreg --pedir="$RunPEDIR" ##> /dev/null 2>&1
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r"_warp.nii.gz

			# Use BBRegister (BBR) to refine the existing coregistration and write an FSL-style transform.
			bbregister --s freesurfer --mov "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --init-reg "$1"/res0urces/eye.dat --surf white.deformed --bold --reg "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".dat --6 --o "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz ##> /dev/null 2>&1 
			tkregister2 --s freesurfer --noedit --reg "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".dat --mov "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r".nii.gz --targ "$2"/anat/T1w/T1w_acpc_dc_restore.nii.gz --fslregout "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".mat ##> /dev/null 2>&1 

			# add BBR step as post warp linear transformation & generate inverse warp;
			convertwarp --warp1="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg_S"$6"_R"$r"_warp.nii.gz --postmat="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".mat --ref="$5" --out="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			mv "$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r".nii.gz "$COREG_QA_DIR"/SBref2acpc_EpiReg+BBR_ScanSpecificFM_S"$6"_R"$r".nii.gz
			
			# warp SBref image into MNI atlas volume space in a single spline warp; can be used for coregistration QA
			convertwarp --ref="$5" --warp1="$2"/func/xfms/"$FuncXfmsDir"/SBref2acpc_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz --warp2="$2"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz --out="$2"/func/xfms/"$FuncXfmsDir"/SBref2nonlin_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz
			applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --ref="$5" --out="$COREG_QA_DIR"/SBref2nonlin_EpiReg+BBR_ScanSpecificFM_S"$6"_R"$r".nii.gz --warp="$2"/func/xfms/"$FuncXfmsDir"/SBref2nonlin_EpiReg+BBR_S"$6"_R"$r"_warp.nii.gz

		fi

        # repeat warps (ACPC, MNI) but this time with the native --> acpc co-registration using an average field map;
        flirt -dof "$4" -in "$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz -ref "$2"/func/xfms/"$FuncXfmsDir"/AvgSBref.nii.gz -out "$COREG_QA_DIR"/SBref2AvgSBref_S"$6"_R"$r".nii.gz -omat "$COREG_QA_DIR"/SBref2AvgSBref_S"$6"_R"$r".mat
        applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --premat="$COREG_QA_DIR"/SBref2AvgSBref_S"$6"_R"$r".mat --warp="$2"/func/xfms/"$FuncXfmsDir"/AvgSBref2acpc_EpiReg+BBR_warp.nii.gz --out="$COREG_QA_DIR"/SBref2acpc_EpiReg+BBR_AvgFM_S"$6"_R"$r".nii.gz --ref="$5"
        applywarp --interp=spline --in="$2"/func/"$FuncDirName"/session_"$6"/run_"$r"/SBref.nii.gz --premat="$COREG_QA_DIR"/SBref2AvgSBref_S"$6"_R"$r".mat --warp="$2"/func/xfms/"$FuncXfmsDir"/AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz --out="$COREG_QA_DIR"/SBref2nonlin_EpiReg+BBR_AvgFM_S"$6"_R"$r".nii.gz --ref="$5"

	done
}

export FuncDirName FuncXfmsDir WMSEG_NII FM_DIR COREG_QA_DIR
export -f func # also coregister individual SBrefs to the target anatomical image
# Expand the full session array; a scalar array reference selects only the first session.
parallel --jobs "$NTHREADS" func ::: "$MEDIR" ::: "$Subdir" ::: "$Subject" ::: "$DOF" ::: "$AtlasTemplate" ::: "${Sessions[@]}" ::: "$AvgPEDIR"

# Write run-level pointer files needed by later stages.
# (brain mask and subcortical mask in functional space)

# T2w anatomicals (whole brain & brain extracted in Atlas Template space)
if [ -f "$Subdir"/anat/T1w/T2w_acpc_dc_restore.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T2w_acpc_dc_restore.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T2w_acpc_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
else
	echo "[INFO] Skipping T2w functional-grid resample; missing $Subdir/anat/T1w/T2w_acpc_dc_restore.nii.gz"
	rm -f "$XfmsDir"/T2w_acpc_func.nii.gz
fi
if [ -f "$Subdir"/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T2w_acpc_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
else
	echo "[INFO] Skipping T2w brain functional-grid resample; missing $Subdir/anat/T1w/T2w_acpc_dc_restore_brain.nii.gz"
	rm -f "$XfmsDir"/T2w_acpc_brain_func.nii.gz
fi

# T1w anatomicals (whole brain & brain extracted in Atlas Template space)
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T1w_acpc_dc_restore.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T1w_acpc_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
flirt -interp nearestneighbour -in "$Subdir"/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T1w_acpc_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$XfmsDir"/T1w_acpc_brain_func.nii.gz -bin "$XfmsDir"/T1w_acpc_brain_func_mask.nii.gz

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
		TmpLhRibbon="$XfmsDir"/.tmp_lh.ribbon.nii.gz
		TmpRhRibbon="$XfmsDir"/.tmp_rh.ribbon.nii.gz
		CorticalRibbonSrc="$XfmsDir"/.tmp_CorticalRibbon_auto.nii.gz
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
flirt -interp nearestneighbour -in "$CorticalRibbonSrc" -ref "$AtlasTemplate" -out "$XfmsDir"/CorticalRibbon_acpc_func_mask.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$XfmsDir"/CorticalRibbon_acpc_func_mask.nii.gz -bin "$XfmsDir"/CorticalRibbon_acpc_func_mask.nii.gz

# MNINonlinear functional-grid ribbon (uses nonlinear ACPC->standard warp).
NonlinWarp="$Subdir"/anat/MNINonLinear/xfms/acpc_dc2standard.nii.gz
if [ -f "$NonlinWarp" ]; then
	applywarp --interp=nn --in="$CorticalRibbonSrc" --ref="$AtlasTemplate" --warp="$NonlinWarp" --out="$XfmsDir"/CorticalRibbon_nonlin_func_mask.nii.gz
	fslmaths "$XfmsDir"/CorticalRibbon_nonlin_func_mask.nii.gz -bin "$XfmsDir"/CorticalRibbon_nonlin_func_mask.nii.gz
elif [[ "$AtlasSpace" == "MNINonlinear" ]]; then
	echo "ERROR: missing nonlinear warp required for AtlasSpace=MNINonlinear: $NonlinWarp"
	exit 2
else
	echo "[INFO] Nonlinear ribbon mask skipped (missing warp: $NonlinWarp)"
fi
rm -f "$XfmsDir"/.tmp_CorticalRibbon_auto.nii.gz

# MNINonlinear anatomicals (whole brain & brain extracted in Atlas Template space)
if [ -f "$Subdir"/anat/MNINonLinear/T2w_restore.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T2w_restore.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T2w_nonlin_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
if [ -f "$Subdir"/anat/MNINonLinear/T2w_restore_brain.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T2w_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T2w_nonlin_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
if [ -f "$Subdir"/anat/MNINonLinear/T1w_restore.nii.gz ]; then
	flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T1w_restore.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T1w_nonlin_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fi
flirt -interp nearestneighbour -in "$Subdir"/anat/MNINonLinear/T1w_restore_brain.nii.gz -ref "$AtlasTemplate" -out "$XfmsDir"/T1w_nonlin_brain_func.nii.gz -applyxfm -init "$MEDIR"/res0urces/ident.mat
fslmaths "$XfmsDir"/T1w_nonlin_brain_func.nii.gz -bin "$XfmsDir"/T1w_nonlin_brain_func_mask.nii.gz

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

PointerLog="$COREG_QA_DIR/CoregPointerSelection.tsv"
echo -e "session\trun\trho_avgfm\trho_scan_specific\tselection" > "$PointerLog"

for s in "${Sessions[@]}" ; do
	mapfile -t runs < <(
		find -L "$Subdir"/func/"$FuncDirName"/session_"$s" -mindepth 1 -maxdepth 1 -type d -name 'run_*' -printf '%f\n' |
			sed 's/^run_//' |
			sort -n
	)
	for r in "${runs[@]}" ; do
		AvgWarp="$Subdir/func/xfms/${FuncXfmsDir}/AvgSBref2${PointerWarpSpace}_EpiReg+BBR_warp.nii.gz"
		ScanWarp="$Subdir/func/xfms/${FuncXfmsDir}/SBref2${PointerWarpSpace}_EpiReg+BBR_S${s}_R${r}_warp.nii.gz"
		RunSBref="$Subdir/func/$FuncDirName/session_${s}/run_${r}/SBref.nii.gz"
		AvgSBref="$Subdir/func/xfms/${FuncXfmsDir}/AvgSBref.nii.gz"
		TargetTxt="$Subdir/func/$FuncDirName/session_${s}/run_${r}/IntermediateCoregTarget.txt"
		WarpTxt="$Subdir/func/$FuncDirName/session_${s}/run_${r}/Intermediate2ACPCWarp.txt"
		MatTxt="$Subdir/func/$FuncDirName/session_${s}/run_${r}/Intermediate2ACPCMat.txt"

		rm -f "$MatTxt"
		if [[ "$ScanSpecificFM" == "1" && -f "$ScanWarp" ]]; then
			echo "$RunSBref" > "$TargetTxt"
			echo "$ScanWarp" > "$WarpTxt"
			echo -e "${s}\t${r}\tnan\tnan\tscan_specific" >> "$PointerLog"
		else
			echo "$AvgSBref" > "$TargetTxt"
			echo "$AvgWarp" > "$WarpTxt"
			if [[ "${FUNC_NOFIELDMAP_MODE:-0}" == "1" ]]; then
				echo "$Subdir/func/xfms/${FuncXfmsDir}/AvgSBref2acpc_EpiReg+BBR_full.mat" > "$MatTxt"
			fi
			echo -e "${s}\t${r}\tnan\tnan\tavgfm" >> "$PointerLog"
		fi
	done
done

# Remove the temporary FreeSurfer folder.
rm -rf "$Subdir"/anat/T1w/freesurfer/ 

# Generate subcortical ROIs in the ACPC/nonlinear functional grid.
"$MEDIR"/lib/make_precise_subcortical_labels.sh "$Subdir" "$AtlasTemplate" "$MEDIR"

echo "[INFO] Coreg module complete. MATLAB coregistration QA post-steps are not used in this revised pipeline."

#!/usr/bin/env bash
# MEDIC-aware slice-time/headmotion/anatomical resampling.
#
# Conservative implementation:
#   1) Use warpkit's wk-apply-warp directly for frame-wise MEDIC unwarping.
#   2) Estimate motion on MEDIC-unwarped echo-averaged data.
#   3) Apply motion + unwarped-SBref-to-anatomical warp in the normal FSL pass.
#
# This intentionally avoids pretending that FSL/warpkit warp-field conventions
# are validated for one-step composition before we have a phantom/known-answer
# test. It preserves warpkit MEDIC behavior exactly for the distortion step.

set -euo pipefail
shopt -s nullglob

if [[ "$#" -lt 7 || "$#" -gt 10 ]]; then
  echo "Usage: mefmri_func_headmotion_medic.sh <MEDIR> <Subject> <StudyFolder> <AtlasTemplate> <DOF> <NTHREADS> <StartSession> [AtlasSpace] [FuncDirName] [FuncFilePrefix]" >&2
  exit 2
fi

MEDIR="$1"
Subject="$2"
StudyFolder="$3"
AtlasTemplate="$4"
DOF="$5"
NTHREADS="$6"
StartSession="$7"
AtlasSpace="${8:-${AtlasSpace:-T1w}}"
FuncDirName="${9:-${FUNC_DIRNAME:-rest}}"
FuncFilePrefix="${10:-${FUNC_FILE_PREFIX:-Rest}}"
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-$FuncDirName}"

Subdir="$StudyFolder/$Subject"
WARPKIT_ENV="${WARPKIT_ENV:-warpkit_env}"
WARPKIT_ACTIVATE_MODE="${WARPKIT_ACTIVATE_MODE:-conda_activate}"
WARPKIT_APPLY_WARP_BIN="${WARPKIT_APPLY_WARP_BIN:-wk-apply-warp}"
SBREF_ECHO_COMBINATION="${SBREF_ECHO_COMBINATION:-mean}"
SBREF_MAX_TE_MS="${SBREF_MAX_TE_MS:-60}"
SBREF_T2SMAP_MASK_THR="${SBREF_T2SMAP_MASK_THR:-1}"
SBREF_T2SMAP_THREADS="${SBREF_T2SMAP_THREADS:-1}"
SBREF_KEEP_INTERMEDIATES="${SBREF_KEEP_INTERMEDIATES:-0}"
TEDANA_ENV="${TEDANA_ENV:-mefmri_env}"
TEDANA_ACTIVATE_MODE="${TEDANA_ACTIVATE_MODE:-conda_activate}"
T2SMAP_BIN="${T2SMAP_BIN:-t2smap}"

case "$AtlasSpace" in
  T1w|MNINonlinear) ;;
  *)
    echo "ERROR: mefmri_func_headmotion_medic.sh invalid AtlasSpace='$AtlasSpace' (expected T1w or MNINonlinear)" >&2
    exit 2
    ;;
esac

die() {
  echo "ERROR: $*" >&2
  exit 2
}

resolve_warpkit_bin() {
  local exe="$1"
  if [[ "$exe" == */* ]]; then
    [[ -x "$exe" ]] || die "Configured warpkit executable is not executable: $exe"
    echo "$exe"
    return 0
  fi
  case "$WARPKIT_ACTIVATE_MODE" in
    direct)
      command -v "$exe" || die "Could not find $exe on PATH with WARPKIT_ACTIVATE_MODE=direct"
      ;;
    conda_activate|conda_run)
      command -v conda >/dev/null 2>&1 || die "conda is required for WARPKIT_ACTIVATE_MODE=$WARPKIT_ACTIVATE_MODE"
      local base candidate
      base="$(conda info --base 2>/dev/null || true)"
      [[ -n "$base" ]] || die "Could not resolve conda base path"
      candidate="$base/envs/$WARPKIT_ENV/bin/$exe"
      [[ -x "$candidate" ]] || die "Missing $exe in conda env '$WARPKIT_ENV': $candidate"
      echo "$candidate"
      ;;
    *)
      die "Invalid WARPKIT_ACTIVATE_MODE=$WARPKIT_ACTIVATE_MODE"
      ;;
  esac
}

resolve_tedana_bin() {
  local exe="$1"
  if [[ "$exe" == */* ]]; then
    [[ -x "$exe" ]] || die "Configured tedana executable is not executable: $exe"
    echo "$exe"
    return 0
  fi
  case "$TEDANA_ACTIVATE_MODE" in
    direct)
      command -v "$exe" || die "Could not find $exe on PATH with TEDANA_ACTIVATE_MODE=direct"
      ;;
    conda_activate|conda_run)
      command -v conda >/dev/null 2>&1 || die "conda is required for TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE"
      local base candidate
      base="$(conda info --base 2>/dev/null || true)"
      [[ -n "$base" ]] || die "Could not resolve conda base path"
      candidate="$base/envs/$TEDANA_ENV/bin/$exe"
      [[ -x "$candidate" ]] || die "Missing $exe in conda env '$TEDANA_ENV': $candidate"
      echo "$candidate"
      ;;
    *)
      die "Invalid TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE"
      ;;
  esac
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

APPLY_WARP_BIN="$(resolve_warpkit_bin "$WARPKIT_APPLY_WARP_BIN")"

echo "[medic-headmotion] Subject=$Subject FuncDirName=$FuncDirName AtlasSpace=$AtlasSpace"
echo "[medic-headmotion] warpkit apply=$APPLY_WARP_BIN"

mapfile -t AllScans < <(find -L "$Subdir/func/$FuncDirName" -mindepth 2 -maxdepth 2 -type d -path '*/session_*/run_*' | sort -V | awk -v start="$StartSession" '
  {
    split($0, parts, "/session_");
    split(parts[2], sr, "/run_");
    if (sr[1] + 0 >= start) {
      print "/session_" sr[1] "/run_" sr[2]
    }
  }')
[[ "${#AllScans[@]}" -gt 0 ]] || die "no processed run directories found under $Subdir/func/$FuncDirName"

func() {
  local MEDIR_IN="$1"
  local AtlasTemplate_IN="$2"
  local Subdir_IN="$3"
  local DOF_IN="$4"
  local ScanRel="$5"

  local ApplyN4Bias="${APPLY_N4_BIAS:-0}"
  local KeepMCF="${HEADMOTION_KEEP_MCF:-0}"
  local ReorientFuncToStd="${FUNC_REORIENT_TO_STD:-1}"
  local ProcDir="$Subdir_IN/func/$FuncDirName$ScanRel"
  local UnprocDir="$Subdir_IN/func/unprocessed/$FuncDirName$ScanRel"
  local BrainMask="$Subdir_IN/func/xfms/$FuncXfmsDir/T1w_acpc_brain_func_mask.nii.gz"
  if [[ "$AtlasSpace" == "MNINonlinear" ]]; then
    BrainMask="$Subdir_IN/func/xfms/$FuncXfmsDir/T1w_nonlin_brain_func_mask.nii.gz"
  fi

  local s r DMap PED UnwarpedSBRef IntermediateCoregTarget Intermediate2ACPCWarp
  s="$(basename "$(dirname "$ProcDir")" | sed 's/^session_//')"
  r="$(basename "$ProcDir" | sed 's/^run_//')"

  [[ -f "$ProcDir/TE.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/TE.txt" >&2; return 1; }
  [[ -f "$ProcDir/TR.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/TR.txt" >&2; return 1; }
  [[ -f "$ProcDir/SliceTiming.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/SliceTiming.txt" >&2; return 1; }
  [[ -f "$ProcDir/IntermediateCoregTarget.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/IntermediateCoregTarget.txt" >&2; return 1; }
  [[ -f "$ProcDir/Intermediate2ACPCWarp.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/Intermediate2ACPCWarp.txt" >&2; return 1; }
  [[ -f "$ProcDir/MEDICDisplacementMap.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/MEDICDisplacementMap.txt" >&2; return 1; }
  [[ -f "$ProcDir/MEDICPhaseEncodingAxis.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/MEDICPhaseEncodingAxis.txt" >&2; return 1; }
  [[ -f "$ProcDir/MEDICUnwarpedSBRef.txt" ]] || { echo "[ERROR] medic-headmotion missing $ProcDir/MEDICUnwarpedSBRef.txt" >&2; return 1; }

  DMap="$(cat "$ProcDir/MEDICDisplacementMap.txt")"
  PED="$(cat "$ProcDir/MEDICPhaseEncodingAxis.txt")"
  UnwarpedSBRef="$(cat "$ProcDir/MEDICUnwarpedSBRef.txt")"
  IntermediateCoregTarget="$(cat "$ProcDir/IntermediateCoregTarget.txt")"
  Intermediate2ACPCWarp="$(cat "$ProcDir/Intermediate2ACPCWarp.txt")"

  [[ -f "$DMap" ]] || { echo "[ERROR] medic-headmotion missing displacement map $DMap" >&2; return 1; }
  [[ -f "$UnwarpedSBRef" ]] || { echo "[ERROR] medic-headmotion missing unwarped SBRef $UnwarpedSBRef" >&2; return 1; }
  [[ -f "$IntermediateCoregTarget" ]] || { echo "[ERROR] medic-headmotion missing target $IntermediateCoregTarget" >&2; return 1; }
  [[ -f "$Intermediate2ACPCWarp" ]] || { echo "[ERROR] medic-headmotion missing warp $Intermediate2ACPCWarp" >&2; return 1; }
  [[ -f "$BrainMask" ]] || { echo "[ERROR] medic-headmotion missing brain mask $BrainMask" >&2; return 1; }

  rm -rf "$ProcDir/MCF" "$ProcDir/Rest_AVG_mcf.mat" "$ProcDir/Rest_AVG_mcf.mat+" "$ProcDir/vols"
  rm -f "$ProcDir"/Rest_AVG_distorted.nii.gz "$ProcDir"/Rest_AVG_medic.nii.gz "$ProcDir"/Rest_AVG_mcf.nii.gz "$ProcDir"/Rest_E1_medic.nii.gz
  rm -f "$ProcDir"/"${FuncFilePrefix}"_E*_distorted.nii.gz "$ProcDir"/"${FuncFilePrefix}"_E*_medic.nii.gz "$ProcDir"/"${FuncFilePrefix}"_E*_acpc.nii.gz
  mkdir -p "$ProcDir/vols"

  reorient_to_std_inplace() {
    local img="$1"
    local tmp="${img%.nii.gz}_reorient_tmp.nii.gz"
    [[ "$ReorientFuncToStd" == "1" ]] || return 0
    fslreorient2std "$img" "$tmp"
    mv -f "$tmp" "$img"
  }

  combine_unwarped_echoes() {
    local method="$1"
    local out="$2"
    local tmp_dir="$3"
    shift 3
    local refs=("$@")

    [[ "${#refs[@]}" -gt 0 ]] || { echo "[ERROR] medic-headmotion no echoes selected for scout combination" >&2; return 1; }

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
        local sq
        fslmaths "${refs[0]}" -sqr "$sum"
        local idx
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
        [[ "${#refs[@]}" -ge 2 ]] || { echo "[ERROR] medic-headmotion $method scout combination requires at least two echoes" >&2; return 1; }
        [[ "${#selected_tes_seconds[@]}" -eq "${#refs[@]}" ]] || { echo "[ERROR] medic-headmotion t2smap TE/reference count mismatch" >&2; return 1; }
        local t2smap_bin t2s_dir mask
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
        [[ -f "$t2s_dir/scout_desc-optcom_bold.nii.gz" ]] || { echo "[ERROR] t2smap did not create scout_desc-optcom_bold.nii.gz" >&2; return 1; }
        cp "$t2s_dir/scout_desc-optcom_bold.nii.gz" "$out"
        ;;
      *)
        echo "[ERROR] unsupported SBREF_ECHO_COMBINATION=$method (use mean|sos|t2s|paid|first|last)" >&2
        return 1
        ;;
    esac
  }

  local te tr n_te raw_echo echo_copy raw_e1 nVols rmVols
  te="$(cat "$ProcDir/TE.txt")"
  tr="$(cat "$ProcDir/TR.txt")"
  n_te=0

  raw_e1=("$UnprocDir"/"$FuncFilePrefix"*_E1.nii.gz)
  [[ "${#raw_e1[@]}" -gt 0 ]] || { echo "[ERROR] medic-headmotion missing first-echo input in $UnprocDir" >&2; return 1; }
  nVols="$(fslnvols "${raw_e1[0]}")"
  rmVols=0
  if [[ -f "$ProcDir/rmVols.txt" ]]; then
    rmVols="$(cat "$ProcDir/rmVols.txt")"
  fi

  local EchoWorkDir="$ProcDir/MEDIC/headmotion_echoes"
  rm -rf "$EchoWorkDir"
  mkdir -p "$EchoWorkDir"
  local selected_refs=()
  local selected_tes_seconds=()
  local scout_sources="$ProcDir/MEDIC/headmotion_scout_echo_sources.tsv"
  {
    echo -e "echo\tte_ms\tincluded\tsource_path\tmedic_unwarped_path"
  } > "$scout_sources"

  for echo_time in $te; do
    n_te=$((n_te + 1))
    raw_echo=("$UnprocDir"/"$FuncFilePrefix"*_E"$n_te".nii.gz)
    [[ "${#raw_echo[@]}" -gt 0 ]] || { echo "[ERROR] medic-headmotion missing raw echo $n_te in $UnprocDir" >&2; return 1; }
    echo_copy="$EchoWorkDir/${FuncFilePrefix}_E${n_te}_distorted.nii.gz"
    cp "${raw_echo[0]}" "$echo_copy"
    reorient_to_std_inplace "$echo_copy"

    "$APPLY_WARP_BIN" \
      --input "$echo_copy" \
      --transform "$DMap" \
      --transform-type map \
      --phase-encoding-axis "$PED" \
      --output "$EchoWorkDir/${FuncFilePrefix}_E${n_te}_medic.nii.gz"
    rm -f "$echo_copy"

    local include=1
    if [[ "$SBREF_MAX_TE_MS" != "none" ]]; then
      include="$(float_lt "$echo_time" "$SBREF_MAX_TE_MS")"
    fi
    if [[ "$include" == "1" ]]; then
      selected_refs+=("$EchoWorkDir/${FuncFilePrefix}_E${n_te}_medic.nii.gz")
      selected_tes_seconds+=("$(ms_to_s "$echo_time")")
      echo -e "${n_te}\t${echo_time}\t1\t${raw_echo[0]}\t$EchoWorkDir/${FuncFilePrefix}_E${n_te}_medic.nii.gz" >> "$scout_sources"
    else
      echo -e "${n_te}\t${echo_time}\t0\t${raw_echo[0]}\t$EchoWorkDir/${FuncFilePrefix}_E${n_te}_medic.nii.gz" >> "$scout_sources"
    fi
  done

  combine_unwarped_echoes "$SBREF_ECHO_COMBINATION" "$ProcDir/Rest_AVG_medic.nii.gz" "$EchoWorkDir" "${selected_refs[@]}"
  cp "$EchoWorkDir/${FuncFilePrefix}_E1_medic.nii.gz" "$ProcDir/Rest_E1_medic.nii.gz"

  fslmaths "$ProcDir/Rest_E1_medic.nii.gz" -Tmean "$ProcDir/Mean.nii.gz"
  rm -f "$ProcDir/Rest_E1_medic.nii.gz"
  if [[ "$ApplyN4Bias" -eq 1 ]]; then
    N4BiasFieldCorrection -d 3 -i "$ProcDir/Mean.nii.gz" -o ["$ProcDir/Mean_Restored.nii.gz","$ProcDir/Bias_field.nii.gz"]
    flirt -in "$ProcDir/Bias_field.nii.gz" -ref "$ProcDir/Mean.nii.gz" -applyxfm -init "$MEDIR_IN/res0urces/ident.mat" -out "$ProcDir/Bias_field.nii.gz" -interp spline
    fslmaths "$ProcDir/Rest_AVG_medic.nii.gz" -div "$ProcDir/Bias_field.nii.gz" "$ProcDir/Rest_AVG_medic.nii.gz"
  fi
  rm -f "$ProcDir"/Mean*.nii.gz "$ProcDir"/Bias*.nii.gz

  if [[ "$rmVols" -gt 0 ]]; then
    fslroi "$ProcDir/Rest_AVG_medic.nii.gz" "$ProcDir/Rest_AVG_medic.nii.gz" "$rmVols" $((nVols - rmVols))
  fi

  mcflirt -dof "$DOF_IN" -stages 3 -plots -in "$ProcDir/Rest_AVG_medic.nii.gz" -r "$UnwarpedSBRef" -out "$ProcDir/MCF"
  rm -f "$ProcDir/MCF.nii.gz"

  slicetimer -i "$ProcDir/Rest_AVG_medic.nii.gz" -o "$ProcDir/Rest_AVG_medic.nii.gz" -r "$tr" --tcustom="$ProcDir/SliceTiming.txt"
  mcflirt -dof "$DOF_IN" -mats -stages 3 -in "$ProcDir/Rest_AVG_medic.nii.gz" -r "$IntermediateCoregTarget" -out "$ProcDir/Rest_AVG_mcf"
  rm -f "$ProcDir/Rest_AVG_distorted.nii.gz" "$ProcDir/Rest_AVG_medic.nii.gz"

  for e in $(seq 1 1 "$n_te"); do
    [[ -f "$EchoWorkDir/${FuncFilePrefix}_E${e}_medic.nii.gz" ]] || { echo "[ERROR] medic-headmotion missing unwarped echo $e in $EchoWorkDir" >&2; return 1; }
    cp "$EchoWorkDir/${FuncFilePrefix}_E${e}_medic.nii.gz" "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz"

    if [[ "$rmVols" -gt 0 ]]; then
      fslroi "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz" "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz" "$rmVols" $((nVols - rmVols))
    fi

    slicetimer -i "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz" --tcustom="$ProcDir/SliceTiming.txt" -r "$tr" -o "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz"
    fslsplit "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz" "$ProcDir/Rest_AVG_mcf.mat/vol_" -t

    local mats images
    mats=("$ProcDir"/Rest_AVG_mcf.mat/MAT_*)
    images=("$ProcDir"/Rest_AVG_mcf.mat/vol_*.nii.gz)
    if [[ "${#images[@]}" -ne "${#mats[@]}" ]]; then
      echo "[ERROR] medic-headmotion applywarp pairing mismatch for $ScanRel echo $e: images=${#images[@]} mats=${#mats[@]}" >&2
      return 1
    fi

    for (( i=0; i<${#images[@]}; i++ )); do
      applywarp --interp=spline --in="${images[$i]}" --premat="${mats[$i]}" --warp="$Intermediate2ACPCWarp" --out="${images[$i]}" --ref="$AtlasTemplate_IN"
    done

    fslmerge -t "$ProcDir/${FuncFilePrefix}_E${e}_acpc.nii.gz" "$ProcDir"/Rest_AVG_mcf.mat/*.nii.gz
    fslmaths "$ProcDir/${FuncFilePrefix}_E${e}_acpc.nii.gz" -mas "$BrainMask" "$ProcDir/${FuncFilePrefix}_E${e}_acpc.nii.gz"
    rm -f "$ProcDir"/Rest_AVG_mcf.mat/*.nii.gz "$ProcDir/${FuncFilePrefix}_E${e}_medic.nii.gz"
  done

  rm -rf "$ProcDir/MCF"
  mv "$ProcDir"/*_mcf*.mat "$ProcDir/MCF"

  if [[ "$ApplyN4Bias" -eq 1 ]]; then
    fslmaths "$ProcDir/${FuncFilePrefix}_E1_acpc.nii.gz" -Tmean "$ProcDir/Mean.nii.gz"
    fslmaths "$ProcDir/Mean.nii.gz" -thr 0 "$ProcDir/Mean.nii.gz"
    N4BiasFieldCorrection -d 3 -i "$ProcDir/Mean.nii.gz" -o ["$ProcDir/Mean_Restored.nii.gz","$ProcDir/Bias_field.nii.gz"]
    flirt -in "$ProcDir/Bias_field.nii.gz" -ref "$ProcDir/Mean.nii.gz" -applyxfm -init "$MEDIR_IN/res0urces/ident.mat" -out "$ProcDir/Bias_field.nii.gz" -interp spline
    for e in $(seq 1 1 "$n_te"); do
      fslmaths "$ProcDir/${FuncFilePrefix}_E${e}_acpc.nii.gz" -div "$ProcDir/Bias_field.nii.gz" "$ProcDir/${FuncFilePrefix}_E${e}_acpc.nii.gz"
    done
  fi

  rm -f "$ProcDir"/Mean*.nii.gz "$ProcDir"/Bias*.nii.gz
  if [[ "$KeepMCF" -eq 0 ]]; then
    rm -rf "$ProcDir/MCF" "$ProcDir/Rest_AVG_mcf.mat" "$ProcDir/Rest_AVG_mcf.mat+"
  fi
  if [[ "$SBREF_KEEP_INTERMEDIATES" == "1" ]]; then
    echo "$EchoWorkDir" > "$ProcDir/MEDIC/headmotion_intermediates_dir.txt"
    echo "$ProcDir/vols" > "$ProcDir/MEDIC/headmotion_scout_intermediates_dir.txt"
  else
    rm -rf "$ProcDir/vols"
    rm -rf "$EchoWorkDir"
  fi
}

export FuncDirName FuncFilePrefix FuncXfmsDir AtlasSpace APPLY_WARP_BIN
export SBREF_ECHO_COMBINATION SBREF_MAX_TE_MS SBREF_T2SMAP_MASK_THR SBREF_T2SMAP_THREADS SBREF_KEEP_INTERMEDIATES
export TEDANA_ENV TEDANA_ACTIVATE_MODE T2SMAP_BIN
export -f func die resolve_tedana_bin float_lt ms_to_s
parallel --jobs "$NTHREADS" func ::: "$MEDIR" ::: "$AtlasTemplate" ::: "$Subdir" ::: "$DOF" ::: "${AllScans[@]}"

echo "[medic-headmotion] complete"

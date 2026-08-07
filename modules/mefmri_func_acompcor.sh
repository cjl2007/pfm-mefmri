#!/bin/bash
# Beta aCompCor denoising module.
#
# Current behavior:
# - standalone module, not wired into mefmri_pipeline.sh yet
# - processes each run into a dedicated aCompCor subdir
# - expects a single 4D input timeseries per run plus motion parameters

set -euo pipefail
IFS=$'\n\t'

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage:
  mefmri_func_acompcor.sh <Subject> <StudyFolder> <MEDIR> [StartSession]

Environment overrides:
  ACOMPCOR_INPUT_BASENAME       Input 4D series per run
  ACOMPCOR_MOTION_BASENAME      Motion parameter file per run
  ACOMPCOR_CENSOR_BASENAME      Optional censor vector file per run
  ACOMPCOR_OUT_SUBDIR           Output subdir under each run (default: aCompCor)
  ACOMPCOR_OUTPUT_BASENAME      Main censored-cleaned output basename
  ACOMPCOR_OUTPUT_FULL_BASENAME Optional full-length cleaned output basename
  ACOMPCOR_T2SMAP_SOURCE_SUBDIR Tedana subdir under each run that may contain T2* map
  ACOMPCOR_T2SMAP_SOURCE_NAME   Tedana T2* map filename to copy into aCompCor dir
  ACOMPCOR_T2SMAP_OUTPUT_NAME   Output T2* map filename inside aCompCor dir
  ACOMPCOR_EXTRAAXIAL_STD_PCT   Extra-axial temporal SD threshold in percent
  ACOMPCOR_EMPTY_COMPARTMENT_ACTION ignore|warn|error (default: warn)
  ACOMPCOR_COMPARTMENT_COND_MAX Per-compartment conditioning threshold
  ACOMPCOR_DESIGN_COND_MAX      Final design conditioning threshold
  ACOMPCOR_BANDPASS_ENABLE      Run FSL -bptf bandpass (default: 1)
  ACOMPCOR_BANDPASS_LOW_HZ      Low cutoff (Hz), e.g. 0.01 (default: 0.01)
  ACOMPCOR_BANDPASS_HIGH_HZ     High cutoff (Hz), e.g. 0.10 (default: 0.10)
  ACOMPCOR_BANDPASS_SUFFIX      Temp suffix used during filtering (default: _bp_tmp)
  ACOMPCOR_BPTF_HP_MAGIC        High-pass sigma factor (default: 2)
  ACOMPCOR_BPTF_LP_MAGIC        Low-pass sigma factor (default: 18)
  ACOMPCOR_PUBLISH_OUTPUT       Copy final cleaned output to run root (default: 0)
  ACOMPCOR_PUBLISHED_BASENAME   Run-root output basename when publishing

Example:
  ACOMPCOR_INPUT_BASENAME=Rest_OCME.nii.gz \
  bash modules/mefmri_func_acompcor.sh Sub001 /path/to/study /path/to/wcm-mepipe 1
EOF
  exit 0
fi

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
MEDIR="${3:?missing MEDIR}"
StartSession="${4:-1}"

Subdir="$StudyFolder/$Subject"
FuncDirName="${FUNC_DIRNAME:-rest}"
FuncFilePrefix="${FUNC_FILE_PREFIX:-Rest}"

: "${ACOMPCOR_PYTHON:=${PIPELINE_PYTHON:-python3}}"
: "${ACOMPCOR_OUT_SUBDIR:=aCompCor}"
: "${ACOMPCOR_INPUT_BASENAME:=${FuncFilePrefix}_OCME.nii.gz}"
: "${ACOMPCOR_OUTPUT_BASENAME:=${FuncFilePrefix}_OCME+aCompCor}"
: "${ACOMPCOR_OUTPUT_FULL_BASENAME:=${FuncFilePrefix}_OCME+aCompCor_full}"
: "${ACOMPCOR_MOTION_BASENAME:=MCF.par}"
: "${ACOMPCOR_CENSOR_BASENAME:=}"
: "${ACOMPCOR_CENSOR_KEEP_VALUE:=1}"
: "${ACOMPCOR_T2SMAP_SOURCE_SUBDIR:=Tedana}"
: "${ACOMPCOR_T2SMAP_SOURCE_NAME:=T2starmap.nii.gz}"
: "${ACOMPCOR_T2SMAP_OUTPUT_NAME:=t2smap.nii.gz}"
: "${ACOMPCOR_EXTRAAXIAL_STD_PCT:=2.5}"
: "${ACOMPCOR_EXTRAAXIAL_DILATE_ITERS:=1}"
: "${ACOMPCOR_WM_ERODE_ITERS:=2}"
: "${ACOMPCOR_EMPTY_COMPARTMENT_ACTION:=warn}"
: "${ACOMPCOR_COMPARTMENT_COND_MAX:=30}"
: "${ACOMPCOR_DESIGN_COND_MAX:=250}"
: "${ACOMPCOR_BANDPASS_ENABLE:=1}"
: "${ACOMPCOR_BANDPASS_LOW_HZ:=0.01}"
: "${ACOMPCOR_BANDPASS_HIGH_HZ:=0.10}"
: "${ACOMPCOR_BANDPASS_SUFFIX:=_bp_tmp}"
: "${ACOMPCOR_BPTF_HP_MAGIC:=2}"
: "${ACOMPCOR_BPTF_LP_MAGIC:=18}"
: "${ACOMPCOR_PUBLISH_OUTPUT:=0}"
: "${ACOMPCOR_PUBLISHED_BASENAME:=}"

PY_SCRIPT="$MEDIR/lib/acompcor_denoise.py"
if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "ERROR: missing aCompCor Python backend: $PY_SCRIPT" >&2
  exit 2
fi

bandpass_with_fsl_bptf() {
  local in_nii="$1"
  local out_nii="$2"
  local tr low_hz high_hz hp_magic lp_magic hp_sigma lp_sigma

  if [[ ! -f "$in_nii" ]]; then
    echo "ERROR: missing bandpass input: $in_nii" >&2
    return 2
  fi
  command -v fslmaths >/dev/null 2>&1 || { echo "ERROR: ACOMPCOR_BANDPASS_ENABLE=1 requires fslmaths on PATH" >&2; return 2; }
  command -v fslval >/dev/null 2>&1 || { echo "ERROR: ACOMPCOR_BANDPASS_ENABLE=1 requires fslval on PATH" >&2; return 2; }

  tr="$(fslval "$in_nii" pixdim4)"
  low_hz="$ACOMPCOR_BANDPASS_LOW_HZ"
  high_hz="$ACOMPCOR_BANDPASS_HIGH_HZ"
  hp_magic="$ACOMPCOR_BPTF_HP_MAGIC"
  lp_magic="$ACOMPCOR_BPTF_LP_MAGIC"

  if [[ -z "$tr" || "$tr" == "0" ]]; then
    echo "ERROR: could not determine TR (pixdim4) from: $in_nii" >&2
    return 2
  fi

  if [[ "$low_hz" == "0" || "$low_hz" == "0.0" ]]; then
    hp_sigma="-1"
  else
    hp_sigma="$(python3 - <<PY
import math
tr=float("${tr}")
f=float("${low_hz}")
magic=float("${hp_magic}")
print("{:.6f}".format(1.0/(magic*f*tr)))
PY
)"
  fi

  if [[ "$high_hz" == "0" || "$high_hz" == "0.0" ]]; then
    lp_sigma="-1"
  else
    lp_sigma="$(python3 - <<PY
import math
tr=float("${tr}")
f=float("${high_hz}")
magic=float("${lp_magic}")
print("{:.6f}".format(1.0/(magic*f*tr)))
PY
)"
  fi

  echo "[acompcor] bandpass fslmaths -bptf hp_sigma=${hp_sigma} lp_sigma=${lp_sigma} (TR=${tr}s)"
  fslmaths "$in_nii" -bptf "$hp_sigma" "$lp_sigma" "$out_nii"
}

echo "[acompcor] subject=$Subject start_session=$StartSession"
echo "[acompcor] func/$FuncDirName input=$ACOMPCOR_INPUT_BASENAME output_dir=$ACOMPCOR_OUT_SUBDIR"

sessions=("$Subdir"/func/"$FuncDirName"/session_*)
sessions=$(seq "$StartSession" 1 "${#sessions[@]}")

for s in $sessions; do
  runs=("$Subdir"/func/"$FuncDirName"/session_"$s"/run_*)
  runs=$(seq 1 1 "${#runs[@]}")

  for r in $runs; do
    run_dir="$Subdir/func/$FuncDirName/session_$s/run_$r"
    input_nii="$run_dir/$ACOMPCOR_INPUT_BASENAME"
    motion_txt="$run_dir/$ACOMPCOR_MOTION_BASENAME"
    out_dir="$run_dir/$ACOMPCOR_OUT_SUBDIR"
    out_base="$out_dir/$ACOMPCOR_OUTPUT_BASENAME"
    out_full_base="$out_dir/$ACOMPCOR_OUTPUT_FULL_BASENAME"
    t2smap_src="$run_dir/$ACOMPCOR_T2SMAP_SOURCE_SUBDIR/$ACOMPCOR_T2SMAP_SOURCE_NAME"
    t2smap_dst="$out_dir/$ACOMPCOR_T2SMAP_OUTPUT_NAME"

    [[ -f "$input_nii" ]] || { echo "ERROR: missing aCompCor input: $input_nii" >&2; exit 2; }
    [[ -f "$motion_txt" ]] || { echo "ERROR: missing aCompCor motion file: $motion_txt" >&2; exit 2; }

    mkdir -p "$out_dir"
    echo "[acompcor] run session_$s/run_$r"

    cmd=(
      "$ACOMPCOR_PYTHON" "$PY_SCRIPT"
      --subdir "$Subdir"
      --input "$input_nii"
      --motion "$motion_txt"
      --output-base "$out_base"
      --output-full-base "$out_full_base"
      --censor-keep-value "$ACOMPCOR_CENSOR_KEEP_VALUE"
      --extraaxial-std-pct "$ACOMPCOR_EXTRAAXIAL_STD_PCT"
      --extraaxial-dilate-iters "$ACOMPCOR_EXTRAAXIAL_DILATE_ITERS"
      --wm-erode-iters "$ACOMPCOR_WM_ERODE_ITERS"
      --empty-compartment-action "$ACOMPCOR_EMPTY_COMPARTMENT_ACTION"
      --compartment-cond-max "$ACOMPCOR_COMPARTMENT_COND_MAX"
      --design-cond-max "$ACOMPCOR_DESIGN_COND_MAX"
    )
    if [[ -n "$ACOMPCOR_CENSOR_BASENAME" ]]; then
      censor_txt="$run_dir/$ACOMPCOR_CENSOR_BASENAME"
      [[ -f "$censor_txt" ]] || { echo "ERROR: missing aCompCor censor file: $censor_txt" >&2; exit 2; }
      cmd+=( --censor "$censor_txt" )
    fi

    "${cmd[@]}"

    if [[ -f "$t2smap_src" ]]; then
      cp -f "$t2smap_src" "$t2smap_dst"
      echo "[acompcor] copied t2smap -> $t2smap_dst"
    else
      echo "[acompcor] t2smap not found, skipping: $t2smap_src"
    fi

    if [[ "$ACOMPCOR_BANDPASS_ENABLE" == "1" ]]; then
      if [[ "$ACOMPCOR_BANDPASS_LOW_HZ" == "0" && "$ACOMPCOR_BANDPASS_HIGH_HZ" == "0" ]]; then
        echo "ERROR: ACOMPCOR_BANDPASS_ENABLE=1 requires ACOMPCOR_BANDPASS_LOW_HZ and/or ACOMPCOR_BANDPASS_HIGH_HZ" >&2
        exit 2
      fi
      bandpass_with_fsl_bptf "${out_base}.nii.gz" "${out_base}${ACOMPCOR_BANDPASS_SUFFIX}.nii.gz"
      mv -f "${out_base}${ACOMPCOR_BANDPASS_SUFFIX}.nii.gz" "${out_base}.nii.gz"
      if [[ -f "${out_full_base}.nii.gz" ]]; then
        bandpass_with_fsl_bptf "${out_full_base}.nii.gz" "${out_full_base}${ACOMPCOR_BANDPASS_SUFFIX}.nii.gz"
        mv -f "${out_full_base}${ACOMPCOR_BANDPASS_SUFFIX}.nii.gz" "${out_full_base}.nii.gz"
      fi
    fi

    if [[ "$ACOMPCOR_PUBLISH_OUTPUT" == "1" ]]; then
      [[ -n "$ACOMPCOR_PUBLISHED_BASENAME" ]] || { echo "ERROR: ACOMPCOR_PUBLISH_OUTPUT=1 requires ACOMPCOR_PUBLISHED_BASENAME" >&2; exit 2; }
      cp -f "${out_base}.nii.gz" "$run_dir/${ACOMPCOR_PUBLISHED_BASENAME}.nii.gz"
      echo "[acompcor] published cleaned output -> $run_dir/${ACOMPCOR_PUBLISHED_BASENAME}.nii.gz"
    fi
  done
done

echo "[acompcor] done"

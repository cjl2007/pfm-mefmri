#!/usr/bin/env bash
# Generate per-run MEDIC field maps/displacement maps with warpkit.
#
# MEDIC requires per-echo phase companions for each functional echo. The
# pipeline runner uses these outputs for the MEDIC-aware coregistration and
# headmotion paths when DISTORTION_CORRECTION_MODE=medic.

set -euo pipefail
shopt -s nullglob

if [[ "$#" -lt 6 || "$#" -gt 7 ]]; then
  echo "Usage: mefmri_func_medic.sh <MEDIR> <Subject> <StudyFolder> <NTHREADS> <StartSession> <FuncDirName> [FuncFilePrefix]" >&2
  exit 2
fi

MEDIR="$1"
Subject="$2"
StudyFolder="$3"
NTHREADS="$4"
StartSession="$5"
FuncDirName="$6"
FuncFilePrefix="${7:-${FUNC_FILE_PREFIX:-Rest}}"

Subdir="$StudyFolder/$Subject"
UnprocRoot="$Subdir/func/unprocessed/$FuncDirName"
FuncRoot="$Subdir/func/$FuncDirName"
MedicRoot="$FuncRoot/MEDIC"

WARPKIT_ENV="${WARPKIT_ENV:-warpkit_env}"
WARPKIT_ACTIVATE_MODE="${WARPKIT_ACTIVATE_MODE:-conda_activate}"
WARPKIT_MEDIC_BIN="${WARPKIT_MEDIC_BIN:-wk-medic}"
WARPKIT_COMPUTE_JACOBIAN_BIN="${WARPKIT_COMPUTE_JACOBIAN_BIN:-wk-compute-jacobian}"
MEDIC_OVERWRITE="${MEDIC_OVERWRITE:-0}"
MEDIC_NOISE_FRAMES="${MEDIC_NOISE_FRAMES:-0}"
MEDIC_DEBUG="${MEDIC_DEBUG:-0}"
MEDIC_WRAP_LIMIT="${MEDIC_WRAP_LIMIT:-0}"
MEDIC_COMPUTE_JACOBIAN="${MEDIC_COMPUTE_JACOBIAN:-1}"

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
      die "Invalid WARPKIT_ACTIVATE_MODE=$WARPKIT_ACTIVATE_MODE (use conda_activate|conda_run|direct)"
      ;;
  esac
}

json_value() {
  local json="$1"
  local key="$2"
  python3 - "$json" "$key" <<'PY'
import json
import sys
with open(sys.argv[1]) as f:
    data = json.load(f)
value = data.get(sys.argv[2], "")
print(value)
PY
}

[[ -d "$UnprocRoot" ]] || die "Missing unprocessed functional root: $UnprocRoot"
[[ "$MEDIC_NOISE_FRAMES" =~ ^[0-9]+$ ]] || die "MEDIC_NOISE_FRAMES must be a non-negative integer"

MEDIC_BIN="$(resolve_warpkit_bin "$WARPKIT_MEDIC_BIN")"
JACOBIAN_BIN=""
if [[ "$MEDIC_COMPUTE_JACOBIAN" == "1" ]]; then
  JACOBIAN_BIN="$(resolve_warpkit_bin "$WARPKIT_COMPUTE_JACOBIAN_BIN")"
fi

mkdir -p "$MedicRoot"
Manifest="$MedicRoot/medic_manifest.tsv"
echo -e "session\trun\tstatus\tmagnitude_count\tphase_count\tmetadata_count\tout_prefix\tmessage" > "$Manifest"

echo "[medic] Subject=$Subject FuncDirName=$FuncDirName FuncFilePrefix=$FuncFilePrefix"
echo "[medic] warpkit env=$WARPKIT_ENV mode=$WARPKIT_ACTIVATE_MODE wk-medic=$MEDIC_BIN"
echo "[medic] manifest root: $MedicRoot"

mapfile -t SessionDirs < <(find -L "$UnprocRoot" -mindepth 1 -maxdepth 1 -type d -name 'session_*' | sort -V)
[[ "${#SessionDirs[@]}" -gt 0 ]] || die "No session_* directories under $UnprocRoot"

for SessionDir in "${SessionDirs[@]}"; do
  s="${SessionDir##*/session_}"
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  (( s >= StartSession )) || continue

  mapfile -t RunDirs < <(find -L "$SessionDir" -mindepth 1 -maxdepth 1 -type d -name 'run_*' | sort -V)
  for RunDir in "${RunDirs[@]}"; do
    r="${RunDir##*/run_}"
    [[ "$r" =~ ^[0-9]+$ ]] || continue

    mapfile -t Mags < <(find -L "$RunDir" -maxdepth 1 -type f -name "${FuncFilePrefix}_S${s}_R${r}_E*.nii.gz" ! -name '*_phase.nii.gz' | sort -V)
    Phases=()
    Jsons=()
    for mag in "${Mags[@]}"; do
      phase="${mag%.nii.gz}_phase.nii.gz"
      json="${mag%.nii.gz}.json"
      [[ -f "$phase" ]] || die "Missing phase companion for $mag: $phase"
      [[ -f "$json" ]] || die "Missing JSON sidecar for $mag: $json"
      Phases+=("$phase")
      Jsons+=("$json")
    done

    if [[ "${#Mags[@]}" -lt 2 ]]; then
      die "MEDIC requires at least two echoes; found ${#Mags[@]} in $RunDir"
    fi

    OutDir="$FuncRoot/session_${s}/run_${r}/MEDIC"
    mkdir -p "$OutDir"
    OutPrefix="$OutDir/${FuncFilePrefix}_S${s}_R${r}"
    Expected="$OutPrefix"_displacementmaps.nii

    if [[ -f "$Expected" && "$MEDIC_OVERWRITE" != "1" ]]; then
      echo "[medic] skip existing session_${s}/run_${r}: $Expected"
      echo -e "${s}\t${r}\tskipped\t${#Mags[@]}\t${#Phases[@]}\t${#Jsons[@]}\t${OutPrefix}\texisting output" >> "$Manifest"
      continue
    fi

    Cmd=(
      "$MEDIC_BIN"
      --magnitude "${Mags[@]}"
      --phase "${Phases[@]}"
      --metadata "${Jsons[@]}"
      --out-prefix "$OutPrefix"
      --n-cpus "$NTHREADS"
    )
    InputsTsv="$OutDir/medic_inputs.tsv"
    echo -e "echo\tmagnitude\tphase\tmetadata" > "$InputsTsv"
    for i in "${!Mags[@]}"; do
      echo -e "$((i + 1))\t${Mags[$i]}\t${Phases[$i]}\t${Jsons[$i]}" >> "$InputsTsv"
    done
    if (( MEDIC_NOISE_FRAMES > 0 )); then
      Cmd+=(--noiseframes "$MEDIC_NOISE_FRAMES")
    fi
    [[ "$MEDIC_DEBUG" == "1" ]] && Cmd+=(--debug)
    [[ "$MEDIC_WRAP_LIMIT" == "1" ]] && Cmd+=(--wrap-limit)

    echo "[medic] run session_${s}/run_${r} echoes=${#Mags[@]}"
    printf '%q ' "${Cmd[@]}" > "$OutDir/wk_medic_command.txt"
    printf '\n' >> "$OutDir/wk_medic_command.txt"
    "${Cmd[@]}"

    [[ -f "$OutPrefix"_fieldmaps_native.nii ]] || die "wk-medic did not create ${OutPrefix}_fieldmaps_native.nii"
    [[ -f "$OutPrefix"_displacementmaps.nii ]] || die "wk-medic did not create ${OutPrefix}_displacementmaps.nii"
    [[ -f "$OutPrefix"_fieldmaps.nii ]] || die "wk-medic did not create ${OutPrefix}_fieldmaps.nii"

    if [[ "$MEDIC_COMPUTE_JACOBIAN" == "1" ]]; then
      PED="$(json_value "${Jsons[0]}" PhaseEncodingDirection)"
      [[ -n "$PED" ]] || die "Missing PhaseEncodingDirection in ${Jsons[0]}"
      "$JACOBIAN_BIN" \
        --input "$OutPrefix"_displacementmaps.nii \
        --from map \
        --axis "$PED" \
        --output "$OutPrefix"_jacobian.nii
    fi

    echo -e "${s}\t${r}\tcreated\t${#Mags[@]}\t${#Phases[@]}\t${#Jsons[@]}\t${OutPrefix}\tok" >> "$Manifest"
  done
done

echo "[medic] complete. Manifest: $Manifest"
echo "[medic] downstream MEDIC coreg/headmotion is handled by the pipeline runner when DISTORTION_CORRECTION_MODE=medic."

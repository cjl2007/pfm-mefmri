#!/usr/bin/env bash
# No-MATLAB fieldmap module.
# Call signature:
#   mefmri_func_fieldmaps.sh <MEDIR> <Subject> <StudyFolder> <NTHREADS> <StartSession> [FuncDirName]

set -euo pipefail
IFS=$'\n\t'

MEDIR="${1:?missing MEDIR}"
Subject="${2:?missing Subject}"
StudyFolder="${3:?missing StudyFolder}"
NTHREADS="${4:?missing NTHREADS}"
StartSession="${5:-1}"
FuncDirName="${6:-${FUNC_DIRNAME:-rest}}"

Subdir="$StudyFolder/$Subject"
[[ -d "$Subdir" ]] || { echo "ERROR: missing subject dir: $Subdir" >&2; exit 2; }

# Optional knobs from wrapper config.
FM_RAW_DIR_REL="${FM_RAW_DIR_REL:-func/unprocessed/${FuncDirName}/field_maps}"
FM_OUT_DIR_REL="${FM_OUT_DIR_REL:-func/${FuncDirName}/field_maps}"
QC_OUTPUT_DIRNAME="${QC_OUTPUT_DIRNAME:-QC}"
FM_QA_DIR_REL="${FM_QA_DIR_REL:-func/${FuncDirName}/${QC_OUTPUT_DIRNAME}}"
TOPUP_CONFIG="${TOPUP_CONFIG:-b02b0.cnf}"
FM_BET_FRAC="${FM_BET_FRAC:-0.35}"
FM_SMOOTH_SIGMA_MM="${FM_SMOOTH_SIGMA_MM:-2}"
USE_WB_SMOOTHING="${USE_WB_SMOOTHING:-1}"
CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-1}"
FIELDMAPS_PYTHON="${FIELDMAPS_PYTHON:-${PIPELINE_PYTHON:-python3}}"
FUNC_NOFIELDMAP_MODE="${FUNC_NOFIELDMAP_MODE:-0}"
FIELDMAP_STRATEGY="${FIELDMAP_STRATEGY:-topup}" # topup|direct_b0
FM_DIRECT_B0_DELTA_TE_MS="${FM_DIRECT_B0_DELTA_TE_MS:-2.65}"
FM_DIRECT_B0_PHASE_PATTERN="${FM_DIRECT_B0_PHASE_PATTERN:-FieldMap_S*_R*.nii.gz,*_phasediff.nii.gz}"
FM_DIRECT_B0_MAG_PATTERN="${FM_DIRECT_B0_MAG_PATTERN:-Mag_S*_R*.nii.gz,*_magnitude1.nii.gz}"

# Phase encoding handling:
# - infer from JSON by default
# - allow user override via BIDS string (i/j/k with optional -)
FM_PE_MODE="${FM_PE_MODE:-infer}"          # infer|config
FM_AP_PE_DIR="${FM_AP_PE_DIR:-}"           # e.g. j-
FM_PA_PE_DIR="${FM_PA_PE_DIR:-}"           # e.g. j
FM_DEFAULT_AP_VEC="${FM_DEFAULT_AP_VEC:-0 -1 0}"
FM_DEFAULT_PA_VEC="${FM_DEFAULT_PA_VEC:-0 1 0}"

RAW_FM_DIR="$Subdir/$FM_RAW_DIR_REL"
FM_DIR="$Subdir/$FM_OUT_DIR_REL"
QA_DIR="$Subdir/$FM_QA_DIR_REL"
ALL_DIR="$FM_DIR/AllFMs"
TOPUP_DIR="$ALL_DIR/topup"
LOG_DIR="$FM_DIR/logs"
ACQ="$FM_DIR/acqparams.txt"
FM_MODE_TXT="$FM_DIR/FIELDMAP_MODE.txt"
FM_QA_TXT="$QA_DIR/FieldMapFallback.txt"

mkdir -p "$FM_DIR" "$ALL_DIR" "$TOPUP_DIR" "$LOG_DIR" "$QA_DIR"

T1="$Subdir/anat/T1w/T1w_acpc_dc_restore.nii.gz"
T1B="$Subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz"
ASEG="$Subdir/anat/T1w/$Subject/mri/aparc+aseg.mgz"
export SUBJECTS_DIR="$Subdir/anat/T1w"

EPI_REG_DOF="$MEDIR/res0urces/epi_reg_dof"
EYE_DAT="$MEDIR/res0urces/eye.dat"

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }
cleanup_legacy_txt_artifacts() {
  # Legacy pipeline runs may leave these text files; remove them automatically.
  find "$Subdir" -type f \( -name "AllFM.txt" -o -name "AllFMs.txt" \) -delete 2>/dev/null || true
}
trap cleanup_legacy_txt_artifacts EXIT
[[ "$FUNC_NOFIELDMAP_MODE" == "0" || "$FUNC_NOFIELDMAP_MODE" == "1" ]] || die "FUNC_NOFIELDMAP_MODE must be 0 or 1"
log "Functional task: ${FuncDirName}"
log "Fieldmap strategy: ${FIELDMAP_STRATEGY}"
log "Raw fieldmap dir: ${RAW_FM_DIR}"
log "Processed fieldmap dir: ${FM_DIR}"

for c in "$FIELDMAPS_PYTHON" fslmaths; do need_cmd "$c"; done
WB_OK=0
if [[ "$FUNC_NOFIELDMAP_MODE" == "0" ]]; then
  case "$FIELDMAP_STRATEGY" in
    topup)
      for c in topup mcflirt fslmerge flirt bet fslnvols convert_xfm parallel; do need_cmd "$c"; done
      ;;
    direct_b0)
      for c in fsl_prepare_fieldmap fslreorient2std fslmaths fslmerge flirt bet convert_xfm parallel; do need_cmd "$c"; done
      ;;
    *)
      die "FIELDMAP_STRATEGY must be topup or direct_b0 (got: $FIELDMAP_STRATEGY)"
      ;;
  esac
  for c in mri_binarize mri_convert bbregister tkregister2; do need_cmd "$c"; done
  if command -v wb_command >/dev/null 2>&1 && [[ "$USE_WB_SMOOTHING" == "1" ]]; then
    WB_OK=1
  fi
fi

[[ -f "$T1" ]] || die "Missing: $T1"
[[ -f "$T1B" ]] || die "Missing: $T1B"

write_no_fieldmap_outputs() {
  local reason="$1"
  rm -rf "$ALL_DIR" "$TOPUP_DIR"
  mkdir -p "$ALL_DIR" "$TOPUP_DIR" "$QA_DIR"
  printf 'none\n' > "$FM_MODE_TXT"
  cat > "$FM_QA_TXT" <<EOF
No-fieldmap functional mode was used for this subject.
Reason: $reason
Susceptibility distortion correction was skipped.
Avg_FM_rads_acpc.nii.gz is a zero-valued placeholder in T1w/ACPC space.
Avg_FM_mag_acpc.nii.gz and Avg_FM_mag_acpc_brain.nii.gz are anatomical compatibility images.
EOF
  fslmaths "$T1B" -mul 0 "$FM_DIR/Avg_FM_rads_acpc.nii.gz" >/dev/null 2>&1
  cp -f "$T1" "$FM_DIR/Avg_FM_mag_acpc.nii.gz"
  cp -f "$T1B" "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz"
  : > "$LOG_DIR/topup_parallel.log"
  log "Fieldmap module complete in no-fieldmap mode."
  exit 0
}

if [[ "$FUNC_NOFIELDMAP_MODE" == "1" ]]; then
  write_no_fieldmap_outputs "FUNC_NOFIELDMAP_MODE=1"
fi

[[ -d "$RAW_FM_DIR" ]] || die "Missing raw fieldmap dir: $RAW_FM_DIR"
[[ -f "$ASEG" ]] || die "Missing: $ASEG"
[[ -x "$EPI_REG_DOF" ]] || die "Missing executable: $EPI_REG_DOF"
[[ -f "$EYE_DAT" ]] || die "Missing: $EYE_DAT"
printf 'measured\n' > "$FM_MODE_TXT"

if [[ "$FIELDMAP_STRATEGY" == "direct_b0" ]]; then
  printf 'direct_b0\n' > "$FM_MODE_TXT"

  mapfile -t DIRECT_PHASE_FILES < <(
    IFS=',' read -r -a pats <<< "$FM_DIRECT_B0_PHASE_PATTERN"
    for pat in "${pats[@]}"; do
      find -L "$RAW_FM_DIR" -maxdepth 1 -type f -name "$pat" -print 2>/dev/null
    done | sort -V
  )
  [[ "${#DIRECT_PHASE_FILES[@]}" -gt 0 ]] || die "No direct B0 phase/phasediff files found in $RAW_FM_DIR"

  direct_tag_from_phase() {
    local base="$1"
    base="$(basename "$base")"
    base="${base%.nii.gz}"
    if [[ "$base" =~ FieldMap_(S[0-9]+_R[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
    elif [[ "$base" =~ (S[0-9]+_R[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
    else
      printf '%s\n' "$base"
    fi
  }

  direct_find_mag() {
    local tag="$1"
    local candidate
    for candidate in \
      "$RAW_FM_DIR/Mag_${tag}.nii.gz" \
      "$RAW_FM_DIR/${tag}_magnitude1.nii.gz" \
      "$RAW_FM_DIR/${tag}_magnitude.nii.gz"; do
      [[ -f "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
    done
    IFS=',' read -r -a pats <<< "$FM_DIRECT_B0_MAG_PATTERN"
    for pat in "${pats[@]}"; do
      while IFS= read -r candidate; do
        [[ "$(basename "$candidate")" == *"$tag"* ]] && { printf '%s\n' "$candidate"; return 0; }
      done < <(find -L "$RAW_FM_DIR" -maxdepth 1 -type f -name "$pat" -print 2>/dev/null | sort -V)
    done
    return 1
  }

  DIRECT_TAGS=()
  DIRECT_PHASE_BY_TAG=()
  DIRECT_MAG_BY_TAG=()
  for phase in "${DIRECT_PHASE_FILES[@]}"; do
    tag="$(direct_tag_from_phase "$phase")"
    ses="${tag#S}"; ses="${ses%%_*}"
    if [[ "$ses" =~ ^[0-9]+$ ]] && [[ "$ses" -lt "$StartSession" ]]; then
      continue
    fi
    if mag="$(direct_find_mag "$tag")"; then
      DIRECT_TAGS+=( "$tag" )
      DIRECT_PHASE_BY_TAG+=( "$phase" )
      DIRECT_MAG_BY_TAG+=( "$mag" )
    else
      log "Skipping direct B0 tag without magnitude image: $tag"
    fi
  done
  [[ "${#DIRECT_TAGS[@]}" -gt 0 ]] || die "No direct B0 phase/magnitude pairs found at/after StartSession=$StartSession"
  log "Direct B0 tags: ${DIRECT_TAGS[*]}"

  mri_binarize --i "$ASEG" --wm --o "$Subdir/anat/T1w/$Subject/mri/white.mgz" >/dev/null 2>&1
  mri_convert -i "$Subdir/anat/T1w/$Subject/mri/white.mgz" \
    -o "$Subdir/anat/T1w/$Subject/mri/white.nii.gz" --like "$T1" >/dev/null 2>&1
  rm -rf "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1 || true
  cp -rf "$Subdir/anat/T1w/$Subject" "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1

  direct_b0_one() {
    local phase="$1"
    local mag="$2"
    local tag="$3"
    local fm_dir="$4"
    local subdir="$5"
    local delta_te_ms="$6"
    local outroot="$fm_dir/AllFMs/direct_b0/$tag"
    mkdir -p "$outroot"

    cp -f "$phase" "$outroot/FieldMap_${tag}.nii.gz"
    cp -f "$mag" "$outroot/FM_mag_${tag}.nii.gz"
    fslreorient2std "$outroot/FieldMap_${tag}.nii.gz" "$outroot/FieldMap_${tag}.nii.gz" >/dev/null 2>&1
    fslreorient2std "$outroot/FM_mag_${tag}.nii.gz" "$outroot/FM_mag_${tag}.nii.gz" >/dev/null 2>&1
    fslmaths "$outroot/FM_mag_${tag}.nii.gz" -Tmean "$outroot/FM_mag_${tag}.nii.gz" >/dev/null 2>&1
    bet "$outroot/FM_mag_${tag}.nii.gz" "$outroot/FM_mag_brain_${tag}.nii.gz" -f 0.35 -R >/dev/null 2>&1

    flirt -in "$outroot/FM_mag_brain_${tag}.nii.gz" \
      -ref "$subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
      -omat "$outroot/fm2acpc_${tag}.mat" -dof 6 >/dev/null 2>&1
    convert_xfm -omat "$outroot/fm2acpc_inv_${tag}.mat" -inverse "$outroot/fm2acpc_${tag}.mat" >/dev/null 2>&1
    flirt -dof 6 -interp trilinear \
      -in "$subdir/anat/T1w/T1w_acpc_brain_mask.nii.gz" \
      -ref "$outroot/FM_mag_brain_${tag}.nii.gz" \
      -out "$outroot/Mask.nii.gz" -applyxfm -init "$outroot/fm2acpc_inv_${tag}.mat" >/dev/null 2>&1
    fslmaths "$outroot/FM_mag_${tag}.nii.gz" -mas "$outroot/Mask.nii.gz" "$outroot/FM_mag_brain_${tag}.nii.gz" >/dev/null 2>&1

    fsl_prepare_fieldmap SIEMENS \
      "$outroot/FieldMap_${tag}.nii.gz" \
      "$outroot/FM_mag_brain_${tag}.nii.gz" \
      "$outroot/FM_rads_${tag}.nii.gz" \
      "$delta_te_ms" >/dev/null 2>&1

    flirt -dof 6 -interp trilinear -in "$outroot/FM_mag_${tag}.nii.gz" \
      -ref "$subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
      -out "$outroot/FM_mag_acpc_${tag}.nii.gz" -applyxfm -init "$outroot/fm2acpc_${tag}.mat" >/dev/null 2>&1
    fslmaths "$outroot/FM_mag_acpc_${tag}.nii.gz" \
      -mas "$subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
      "$outroot/FM_mag_acpc_brain_${tag}.nii.gz" >/dev/null 2>&1
    flirt -dof 6 -interp trilinear -in "$outroot/FM_rads_${tag}.nii.gz" \
      -ref "$subdir/anat/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
      -out "$outroot/FM_rads_acpc_${tag}.nii.gz" -applyxfm -init "$outroot/fm2acpc_${tag}.mat" >/dev/null 2>&1
  }
  export -f direct_b0_one

  parallel --jobs "$NTHREADS" direct_b0_one ::: "${DIRECT_PHASE_BY_TAG[@]}" :::+ "${DIRECT_MAG_BY_TAG[@]}" :::+ "${DIRECT_TAGS[@]}" ::: "$FM_DIR" ::: "$Subdir" ::: "$FM_DIRECT_B0_DELTA_TE_MS" \
    >"$LOG_DIR/direct_b0_parallel.log" 2>&1 || die "Direct B0 fieldmap processing failed. See $LOG_DIR/direct_b0_parallel.log"

  for tag in "${DIRECT_TAGS[@]}"; do
    OUTROOT="$ALL_DIR/direct_b0/$tag"
    if [[ "$WB_OK" == "1" ]]; then
      wb_command -volume-smoothing "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" "$FM_SMOOTH_SIGMA_MM" \
        "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" -fix-zeros >/dev/null 2>&1
    fi
    ses="${tag#S}"; ses="${ses%%_*}"
    run="${tag#*_R}"
    mv -f "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" "$ALL_DIR/FM_rads_acpc_S${ses}_R${run}.nii.gz"
    mv -f "$OUTROOT/FM_mag_acpc_${tag}.nii.gz" "$ALL_DIR/FM_mag_acpc_S${ses}_R${run}.nii.gz"
    mv -f "$OUTROOT/FM_mag_acpc_brain_${tag}.nii.gz" "$ALL_DIR/FM_mag_acpc_brain_S${ses}_R${run}.nii.gz"
    [[ "$CLEAN_INTERMEDIATE" == "1" ]] && rm -rf "$OUTROOT"
  done

  mapfile -t RADS < <(find "$ALL_DIR" -maxdepth 1 -type f -name 'FM_rads_acpc_S*_R*.nii.gz' | sort -V)
  mapfile -t MAG < <(find "$ALL_DIR" -maxdepth 1 -type f -name 'FM_mag_acpc_S*_R*.nii.gz' | sort -V)
  mapfile -t MAGB < <(find "$ALL_DIR" -maxdepth 1 -type f -name 'FM_mag_acpc_brain_S*_R*.nii.gz' | sort -V)
  [[ "${#RADS[@]}" -ge 1 && "${#MAG[@]}" -ge 1 ]] || die "No direct B0 scan-specific FM outputs generated"
  fslmerge -t "$FM_DIR/Avg_FM_rads_acpc.nii.gz" "${RADS[@]}" >/dev/null 2>&1
  fslmaths "$FM_DIR/Avg_FM_rads_acpc.nii.gz" -Tmean "$FM_DIR/Avg_FM_rads_acpc.nii.gz" >/dev/null 2>&1
  fslmerge -t "$FM_DIR/Avg_FM_mag_acpc.nii.gz" "${MAG[@]}" >/dev/null 2>&1
  fslmaths "$FM_DIR/Avg_FM_mag_acpc.nii.gz" -Tmean "$FM_DIR/Avg_FM_mag_acpc.nii.gz" >/dev/null 2>&1
  fslmerge -t "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" "${MAGB[@]}" >/dev/null 2>&1
  fslmaths "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" -Tmean "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" >/dev/null 2>&1
  rm -rf "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1 || true
  log "Direct B0 fieldmap module complete."
  exit 0
fi

# Build acqparams.txt and the QA summary from BIDS JSON sidecars.
"$FIELDMAPS_PYTHON" - "$RAW_FM_DIR" "$ACQ" "$QA_DIR/AvgFieldMap.txt" "$FM_PE_MODE" \
  "$FM_AP_PE_DIR" "$FM_PA_PE_DIR" "$FM_DEFAULT_AP_VEC" "$FM_DEFAULT_PA_VEC" <<'PY'
import collections
import json
import math
import sys
from pathlib import Path

raw = Path(sys.argv[1])
acq_out = Path(sys.argv[2])
qa_out = Path(sys.argv[3])
pe_mode = sys.argv[4].strip().lower()
ap_override = sys.argv[5].strip()
pa_override = sys.argv[6].strip()
default_ap_vec = sys.argv[7].strip()
default_pa_vec = sys.argv[8].strip()

def die(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def parse_vec(v):
    parts = v.split()
    if len(parts) != 3:
        die(f"Invalid PE vector '{v}'")
    return [int(float(x)) for x in parts]

def bids_dir_to_vec(d):
    if d == "i": return [1, 0, 0]
    if d == "i-": return [-1, 0, 0]
    if d == "j": return [0, 1, 0]
    if d == "j-": return [0, -1, 0]
    if d == "k": return [0, 0, 1]
    if d == "k-": return [0, 0, -1]
    return None

def load_json(p):
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}

def readout(j):
    if "TotalReadoutTime" in j:
        return float(j["TotalReadoutTime"])
    ees = j.get("EffectiveEchoSpacing")
    rpe = j.get("ReconMatrixPE")
    if ees is not None and rpe is not None:
        return float(ees) * (int(rpe) - 1)
    return math.nan

def mean_valid(values):
    vv = [x for x in values if x == x]
    return float(sum(vv) / len(vv)) if vv else math.nan

ap_niis = sorted(raw.glob("AP_S*_R*.nii.gz"))
pa_niis = sorted(raw.glob("PA_S*_R*.nii.gz"))
if not ap_niis or not pa_niis:
    die(f"Missing AP/PA nii.gz in {raw}")

ap_jsons = [raw / (p.name[:-7] + ".json") for p in ap_niis]
pa_jsons = [raw / (p.name[:-7] + ".json") for p in pa_niis]

ap_ro, pa_ro = [], []
ap_dirs, pa_dirs = [], []
lines = [
    f"Number of AP Field Maps: {len(ap_niis)}",
    f"Number of PA Field Maps: {len(pa_niis)}",
]

for apn in ap_niis:
    tag = apn.name[:-7]
    pan = raw / ("PA_" + tag[3:] + ".nii.gz")
    if not pan.exists():
        lines.append(f"Pair ?: {tag} [MISSING_PA] + PA_{tag[3:]} [MISSING]")
        continue
    apj = load_json(raw / (tag + ".json"))
    paj = load_json(raw / ("PA_" + tag[3:] + ".json"))
    apd = str(apj.get("PhaseEncodingDirection", "Unknown"))
    pad = str(paj.get("PhaseEncodingDirection", "Unknown"))
    ap_dirs.append(apd)
    pa_dirs.append(pad)
    ap_ro.append(readout(apj))
    pa_ro.append(readout(paj))
    lines.append(f"Pair {tag[2:]} : {tag} [{apd}] + PA_{tag[3:]} [{pad}]")

ap_mean = mean_valid(ap_ro)
pa_mean = mean_valid(pa_ro)

def choose_dir(dirs, label):
    c = collections.Counter([d for d in dirs if d in {"i","i-","j","j-","k","k-"}])
    if not c:
        return None
    return c.most_common(1)[0][0]

if pe_mode == "config":
    ap_dir = ap_override or None
    pa_dir = pa_override or None
else:
    ap_dir = ap_override or choose_dir(ap_dirs, "AP")
    pa_dir = pa_override or choose_dir(pa_dirs, "PA")

ap_vec = bids_dir_to_vec(ap_dir) if ap_dir else parse_vec(default_ap_vec)
pa_vec = bids_dir_to_vec(pa_dir) if pa_dir else parse_vec(default_pa_vec)

if ap_mean != ap_mean and pa_mean != pa_mean:
    die("Could not infer TotalReadoutTime from JSONs; add TotalReadoutTime or EffectiveEchoSpacing+ReconMatrixPE")
if ap_mean != ap_mean:
    ap_mean = pa_mean
if pa_mean != pa_mean:
    pa_mean = ap_mean

acq_out.parent.mkdir(parents=True, exist_ok=True)
acq_out.write_text(
    f"{ap_vec[0]} {ap_vec[1]} {ap_vec[2]} {ap_mean}\n"
    f"{pa_vec[0]} {pa_vec[1]} {pa_vec[2]} {pa_mean}\n"
)
lines.append(f"acqparams AP row: {ap_vec[0]} {ap_vec[1]} {ap_vec[2]} {ap_mean}")
lines.append(f"acqparams PA row: {pa_vec[0]} {pa_vec[1]} {pa_vec[2]} {pa_mean}")
qa_out.parent.mkdir(parents=True, exist_ok=True)
qa_out.write_text("\n".join(lines) + "\n")
print(f"Wrote {acq_out}")
print(f"Wrote {qa_out}")
PY

mapfile -t AP_FILES < <(ls -1 "$RAW_FM_DIR"/AP_S*_R*.nii.gz 2>/dev/null | sort || true)
TAGS=()
for ap in "${AP_FILES[@]}"; do
  tag="$(basename "$ap")"
  tag="${tag#AP_}"
  tag="${tag%.nii.gz}"
  pa="$RAW_FM_DIR/PA_${tag}.nii.gz"
  [[ -f "$pa" ]] || continue
  ses="${tag#S}"; ses="${ses%%_*}"
  if [[ "$ses" -ge "$StartSession" ]]; then
    TAGS+=( "$tag" )
  fi
done
[[ "${#TAGS[@]}" -gt 0 ]] || die "No AP/PA tags found at/after StartSession=$StartSession"

log "Fieldmap tags: ${TAGS[*]}"

# WM seg + temp freesurfer alias
mri_binarize --i "$ASEG" --wm --o "$Subdir/anat/T1w/$Subject/mri/white.mgz" >/dev/null 2>&1
mri_convert -i "$Subdir/anat/T1w/$Subject/mri/white.mgz" \
  -o "$Subdir/anat/T1w/$Subject/mri/white.nii.gz" --like "$T1" >/dev/null 2>&1
rm -rf "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1 || true
cp -rf "$Subdir/anat/T1w/$Subject" "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1

topup_one() {
  local raw_fm_dir="$1"
  local fm_dir="$2"
  local acq="$3"
  local tag="$4"
  local TOPUP_CONFIG="$5"
  local OUTROOT="$fm_dir/AllFMs/topup/$tag"
  mkdir -p "$OUTROOT"

  cp -f "$raw_fm_dir/AP_${tag}.nii.gz" "$OUTROOT/AP_${tag}.nii.gz"
  cp -f "$raw_fm_dir/PA_${tag}.nii.gz" "$OUTROOT/PA_${tag}.nii.gz"

  local nVols
  nVols="$(fslnvols "$OUTROOT/AP_${tag}.nii.gz")"
  if [[ "$nVols" -gt 1 ]]; then
    mcflirt -in "$OUTROOT/AP_${tag}.nii.gz" -out "$OUTROOT/AP_${tag}.nii.gz" >/dev/null 2>&1
    fslmaths "$OUTROOT/AP_${tag}.nii.gz" -Tmean "$OUTROOT/AP_${tag}.nii.gz" >/dev/null 2>&1
    mcflirt -in "$OUTROOT/PA_${tag}.nii.gz" -out "$OUTROOT/PA_${tag}.nii.gz" >/dev/null 2>&1
    fslmaths "$OUTROOT/PA_${tag}.nii.gz" -Tmean "$OUTROOT/PA_${tag}.nii.gz" >/dev/null 2>&1
  fi

  fslmerge -t "$OUTROOT/AP_PA_${tag}.nii.gz" "$OUTROOT/AP_${tag}.nii.gz" "$OUTROOT/PA_${tag}.nii.gz" >/dev/null 2>&1
  topup --imain="$OUTROOT/AP_PA_${tag}.nii.gz" \
    --datain="$acq" \
    --iout="$OUTROOT/FM_mag_${tag}.nii.gz" \
    --fout="$OUTROOT/FM_hz_${tag}.nii.gz" \
    --config="$TOPUP_CONFIG" >/dev/null 2>&1
  fslmaths "$OUTROOT/FM_hz_${tag}.nii.gz" -mul 6.283 "$OUTROOT/FM_rads_${tag}.nii.gz" >/dev/null 2>&1
  fslmaths "$OUTROOT/FM_mag_${tag}.nii.gz" -Tmean "$OUTROOT/FM_mag_${tag}.nii.gz" >/dev/null 2>&1
}
export -f topup_one

parallel --jobs "$NTHREADS" topup_one ::: "$RAW_FM_DIR" ::: "$FM_DIR" ::: "$ACQ" ::: "${TAGS[@]}" ::: "$TOPUP_CONFIG" \
  >"$LOG_DIR/topup_parallel.log" 2>&1 || die "TOPUP failed. See $LOG_DIR/topup_parallel.log"

for tag in "${TAGS[@]}"; do
  OUTROOT="$TOPUP_DIR/$tag"
  bet "$OUTROOT/FM_mag_${tag}.nii.gz" "$OUTROOT/FM_mag_brain_${tag}.nii.gz" -f "$FM_BET_FRAC" -R >/dev/null 2>&1

  "$EPI_REG_DOF" --epi="$OUTROOT/FM_mag_${tag}.nii.gz" \
    --t1="$T1" --t1brain="$T1B" --out="$OUTROOT/fm2acpc_${tag}" \
    --wmseg="$Subdir/anat/T1w/$Subject/mri/white.nii.gz" --dof=6 >/dev/null 2>&1

  bbregister --s freesurfer --mov "$OUTROOT/fm2acpc_${tag}.nii.gz" \
    --init-reg "$EYE_DAT" --surf white.deformed --bold \
    --reg "$OUTROOT/fm2acpc_bbr_${tag}.dat" --6 \
    --o "$OUTROOT/fm2acpc_bbr_${tag}.nii.gz" >/dev/null 2>&1
  tkregister2 --s freesurfer --noedit --reg "$OUTROOT/fm2acpc_bbr_${tag}.dat" \
    --mov "$OUTROOT/fm2acpc_${tag}.nii.gz" --targ "$T1" \
    --fslregout "$OUTROOT/fm2acpc_bbr_${tag}.mat" >/dev/null 2>&1
  convert_xfm -omat "$OUTROOT/fm2acpc_${tag}.mat" \
    -concat "$OUTROOT/fm2acpc_bbr_${tag}.mat" "$OUTROOT/fm2acpc_${tag}.mat" >/dev/null 2>&1

  flirt -dof 6 -interp spline -in "$OUTROOT/FM_mag_${tag}.nii.gz" -ref "$T1B" \
    -out "$OUTROOT/FM_mag_acpc_${tag}.nii.gz" -applyxfm -init "$OUTROOT/fm2acpc_${tag}.mat" >/dev/null 2>&1
  fslmaths "$OUTROOT/FM_mag_acpc_${tag}.nii.gz" -mas "$T1B" "$OUTROOT/FM_mag_acpc_brain_${tag}.nii.gz" >/dev/null 2>&1
  flirt -dof 6 -interp spline -in "$OUTROOT/FM_rads_${tag}.nii.gz" -ref "$T1B" \
    -out "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" -applyxfm -init "$OUTROOT/fm2acpc_${tag}.mat" >/dev/null 2>&1

  if [[ "$WB_OK" == "1" ]]; then
    wb_command -volume-smoothing "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" "$FM_SMOOTH_SIGMA_MM" \
      "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" -fix-zeros >/dev/null 2>&1
  fi

  ses="${tag#S}"; ses="${ses%%_*}"
  run="${tag#*_R}"
  mv -f "$OUTROOT/FM_rads_acpc_${tag}.nii.gz" "$ALL_DIR/FM_rads_acpc_S${ses}_R${run}.nii.gz"
  mv -f "$OUTROOT/FM_mag_acpc_${tag}.nii.gz" "$ALL_DIR/FM_mag_acpc_S${ses}_R${run}.nii.gz"
  mv -f "$OUTROOT/FM_mag_acpc_brain_${tag}.nii.gz" "$ALL_DIR/FM_mag_acpc_brain_S${ses}_R${run}.nii.gz"
  if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
    rm -rf "$OUTROOT"
  fi
done

mapfile -t RADS < <(ls -1 "$ALL_DIR"/FM_rads_acpc_S*_R*.nii.gz 2>/dev/null | sort || true)
mapfile -t MAG < <(ls -1 "$ALL_DIR"/FM_mag_acpc_S*_R*.nii.gz 2>/dev/null | sort || true)
mapfile -t MAGB < <(ls -1 "$ALL_DIR"/FM_mag_acpc_brain_S*_R*.nii.gz 2>/dev/null | sort || true)
[[ "${#RADS[@]}" -ge 1 && "${#MAG[@]}" -ge 1 ]] || die "No scan-specific FM outputs generated"

fslmerge -t "$FM_DIR/Avg_FM_rads_acpc.nii.gz" "${RADS[@]}" >/dev/null 2>&1
fslmaths "$FM_DIR/Avg_FM_rads_acpc.nii.gz" -Tmean "$FM_DIR/Avg_FM_rads_acpc.nii.gz" >/dev/null 2>&1
fslmerge -t "$FM_DIR/Avg_FM_mag_acpc.nii.gz" "${MAG[@]}" >/dev/null 2>&1
fslmaths "$FM_DIR/Avg_FM_mag_acpc.nii.gz" -Tmean "$FM_DIR/Avg_FM_mag_acpc.nii.gz" >/dev/null 2>&1

if [[ "${#MAGB[@]}" -ge 1 ]]; then
  fslmerge -t "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" "${MAGB[@]}" >/dev/null 2>&1
  fslmaths "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" -Tmean "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" >/dev/null 2>&1
else
  fslmaths "$FM_DIR/Avg_FM_mag_acpc.nii.gz" -mas "$T1B" "$FM_DIR/Avg_FM_mag_acpc_brain.nii.gz" >/dev/null 2>&1
fi

rm -rf "$Subdir/anat/T1w/freesurfer" >/dev/null 2>&1 || true
log "Fieldmap module complete."

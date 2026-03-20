#!/usr/bin/env bash
# No-MATLAB fieldmap module.
# Call signature:
#   mefmri_func_fieldmaps.sh <MEDIR> <Subject> <StudyFolder> <NTHREADS> <StartSession>

set -euo pipefail
IFS=$'\n\t'

MEDIR="${1:?missing MEDIR}"
Subject="${2:?missing Subject}"
StudyFolder="${3:?missing StudyFolder}"
NTHREADS="${4:?missing NTHREADS}"
StartSession="${5:-1}"

Subdir="$StudyFolder/$Subject"
[[ -d "$Subdir" ]] || { echo "ERROR: missing subject dir: $Subdir" >&2; exit 2; }

# Optional knobs from wrapper config.
FM_RAW_DIR_REL="${FM_RAW_DIR_REL:-func/unprocessed/field_maps}"
FM_OUT_DIR_REL="${FM_OUT_DIR_REL:-func/field_maps}"
FM_QA_DIR_REL="${FM_QA_DIR_REL:-func/qa}"
TOPUP_CONFIG="${TOPUP_CONFIG:-b02b0.cnf}"
FM_BET_FRAC="${FM_BET_FRAC:-0.35}"
FM_SMOOTH_SIGMA_MM="${FM_SMOOTH_SIGMA_MM:-2}"
USE_WB_SMOOTHING="${USE_WB_SMOOTHING:-1}"
CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-1}"

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

for c in python3 topup mcflirt fslmaths fslmerge flirt bet fslnvols convert_xfm parallel; do need_cmd "$c"; done
for c in mri_binarize mri_convert bbregister tkregister2; do need_cmd "$c"; done
WB_OK=0
if command -v wb_command >/dev/null 2>&1 && [[ "$USE_WB_SMOOTHING" == "1" ]]; then
  WB_OK=1
fi

[[ -d "$RAW_FM_DIR" ]] || die "Missing raw fieldmap dir: $RAW_FM_DIR"
[[ -f "$T1" ]] || die "Missing: $T1"
[[ -f "$T1B" ]] || die "Missing: $T1B"
[[ -f "$ASEG" ]] || die "Missing: $ASEG"
[[ -x "$EPI_REG_DOF" ]] || die "Missing executable: $EPI_REG_DOF"
[[ -f "$EYE_DAT" ]] || die "Missing: $EYE_DAT"

# Build acqparams.txt and the QA summary from BIDS JSON sidecars.
python3 - "$RAW_FM_DIR" "$ACQ" "$QA_DIR/AvgFieldMap.txt" "$FM_PE_MODE" \
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
  local Subdir="$1"
  local tag="$2"
  local TOPUP_CONFIG="$3"
  local OUTROOT="$Subdir/func/field_maps/AllFMs/topup/$tag"
  mkdir -p "$OUTROOT"

  cp -f "$Subdir/func/unprocessed/field_maps/AP_${tag}.nii.gz" "$OUTROOT/AP_${tag}.nii.gz"
  cp -f "$Subdir/func/unprocessed/field_maps/PA_${tag}.nii.gz" "$OUTROOT/PA_${tag}.nii.gz"

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
    --datain="$Subdir/func/field_maps/acqparams.txt" \
    --iout="$OUTROOT/FM_mag_${tag}.nii.gz" \
    --fout="$OUTROOT/FM_hz_${tag}.nii.gz" \
    --config="$TOPUP_CONFIG" >/dev/null 2>&1
  fslmaths "$OUTROOT/FM_hz_${tag}.nii.gz" -mul 6.283 "$OUTROOT/FM_rads_${tag}.nii.gz" >/dev/null 2>&1
  fslmaths "$OUTROOT/FM_mag_${tag}.nii.gz" -Tmean "$OUTROOT/FM_mag_${tag}.nii.gz" >/dev/null 2>&1
}
export -f topup_one

parallel --jobs "$NTHREADS" topup_one ::: "$Subdir" ::: "${TAGS[@]}" ::: "$TOPUP_CONFIG" \
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

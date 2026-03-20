#!/usr/bin/env bash
# CJL; (cjl2007@med.cornell.edu)
# Tedana-first ME-ICA module with Python-based post-processing.

set -euo pipefail
IFS=$'\n\t'

WORKER_MODE=0
if [[ "${1:-}" == "--worker" ]]; then
  WORKER_MODE=1
  shift
fi

Subject="${1:?missing Subject}"
StudyFolder="${2:?missing StudyFolder}"
NTHREADS="${3:?missing NTHREADS}"
MEPCA="${4:?missing MEPCA}"
MaxIterations="${5:?missing MaxIterations}"
MaxRestarts="${6:?missing MaxRestarts}"
StartSession="${7:?missing StartSession}"
MEDIR="${8:?missing MEDIR}"
WORKER_RUNDIR="${9:-}"

Subdir="$StudyFolder/$Subject"
FuncDirName="${FUNC_DIRNAME:-rest}"
FuncFilePrefix="${FUNC_FILE_PREFIX:-Rest}"
# Coreg currently writes transform-space derivatives under func/xfms/rest.
FuncXfmsDir="${FUNC_XFMS_DIRNAME:-rest}"

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

# Tedana runtime config (wrapper config may override these env vars)
TEDANA_ENV="${TEDANA_ENV:-mefmri_env}"
TEDANA_ACTIVATE_MODE="${TEDANA_ACTIVATE_MODE:-conda_activate}" # conda_activate|conda_run|direct
TEDANA_COMPAT_MODE="${TEDANA_COMPAT_MODE:-modern}"             # modern|auto|legacy|v12|v10
TEDANA_FITTYPE="${TEDANA_FITTYPE:-curvefit}"                   # curvefit|loglin
TEDANA_ICA_METHOD="${TEDANA_ICA_METHOD:-fastica}"              # fastica|robustica
TEDANA_N_ROBUST_RUNS="${TEDANA_N_ROBUST_RUNS:-10}"
TEDANA_SEED="${TEDANA_SEED:-42}"
TEDANA_THREADS="${TEDANA_THREADS:-$NTHREADS}"
TEDANA_MASKTYPE="${TEDANA_MASKTYPE:-none}"                     # none|dropout|decay
TEDANA_CONVENTION="${TEDANA_CONVENTION:-orig}"                 # orig|bids
TEDANA_OVERWRITE="${TEDANA_OVERWRITE:-1}"                      # 0|1
TEDANA_LOWMEM="${TEDANA_LOWMEM:-0}"                            # 0|1
TEDANA_USE_EXTERNAL_MIX="${TEDANA_USE_EXTERNAL_MIX:-0}"        # 0|1
TEDANA_EXTERNAL_MIX_BASENAME="${TEDANA_EXTERNAL_MIX_BASENAME:-}" # e.g., v10_ica_mixing.tsv (looked up in each run dir)
MEICA_PARALLEL_JOBS="${MEICA_PARALLEL_JOBS:-$NTHREADS}"

# Post-tedana reclassification implemented in Python.
MEICA_RECLASSIFY_ENABLE="${MEICA_RECLASSIFY_ENABLE:-1}"       # 0|1
MEICA_CLASSIFIER_MODE="${MEICA_CLASSIFIER_MODE:-nsi}"         # nsi|legacy_template_rho|none
NETWORK_PRIORS_MAT="${NETWORK_PRIORS_MAT:-$MEDIR/res0urces/Priors.mat}"
MEICA_PRIORS_MAT="${MEICA_PRIORS_MAT:-$NETWORK_PRIORS_MAT}"
MEICA_BETAS_CIFTI="${MEICA_BETAS_CIFTI:-}"                    # optional explicit betas_OC.dtseries.nii
MEICA_RHO_RESCUE="${MEICA_RHO_RESCUE:-0.30}"
MEICA_RHO_REJECT="${MEICA_RHO_REJECT:-0.10}"
MEICA_RECLASS_NO_REPORTS="${MEICA_RECLASS_NO_REPORTS:-0}"     # 0|1
MEICA_RECLASS_SUBDIR="${MEICA_RECLASS_SUBDIR:-ManualComponentClassification}"
MEICA_TEDANA_SUBDIR="${MEICA_TEDANA_SUBDIR:-Tedana}"
MEICA_RECLASSIFY_PY="$MEDIR/lib/meica_reclassify_components.py"
MEICA_QC_CIFTI_ENABLE="${MEICA_QC_CIFTI_ENABLE:-1}"           # 0|1
MEICA_QC_CIFTI_TAGS="${MEICA_QC_CIFTI_TAGS:-betas_OC}"        # comma-separated: betas_OC,t2sv,s0v
MEICA_ORIG_ALIAS_ENABLE="${MEICA_ORIG_ALIAS_ENABLE:-1}"       # 1 => expose modern tedana outputs under legacy/orig names
MEICA_NSI_RESCUE_THRESHOLD="${MEICA_NSI_RESCUE_THRESHOLD:-0.20}"
MEICA_NSI_RESCUE_QUANTILE="${MEICA_NSI_RESCUE_QUANTILE:-0.10}"
MEICA_NSI_KILL_MODE="${MEICA_NSI_KILL_MODE:-adaptive}"        # adaptive|fixed
MEICA_NSI_KILL_THRESHOLD="${MEICA_NSI_KILL_THRESHOLD:-0.05}"
MEICA_NSI_KILL_MIN="${MEICA_NSI_KILL_MIN:-0.02}"
MEICA_NSI_KILL_MAX="${MEICA_NSI_KILL_MAX:-0.10}"
MEICA_NSI_KILL_INTERCEPT="${MEICA_NSI_KILL_INTERCEPT:-0.10}"
MEICA_NSI_KILL_SLOPE="${MEICA_NSI_KILL_SLOPE:-0.50}"
MEICA_NSI_GUARDRAIL_KAPPA_RHO="${MEICA_NSI_GUARDRAIL_KAPPA_RHO:-1}" # 0|1
MEICA_SUBCORT_RATIO_THRESH="${MEICA_SUBCORT_RATIO_THRESH:-5.0}"
MEICA_KILL_PRIORITY_ENABLE="${MEICA_KILL_PRIORITY_ENABLE:-1}"       # 0|1
MEICA_KILL_PRIORITY_W_LOGRATIO="${MEICA_KILL_PRIORITY_W_LOGRATIO:-0.50}"
MEICA_KILL_PRIORITY_W_NSI="${MEICA_KILL_PRIORITY_W_NSI:-0.30}"
MEICA_KILL_PRIORITY_W_VAR="${MEICA_KILL_PRIORITY_W_VAR:-0.20}"
MEICA_KILL_VAR_FLOOR_QUANTILE="${MEICA_KILL_VAR_FLOOR_QUANTILE:-0.60}"
MEICA_KILL_CUMVAR_CAP="${MEICA_KILL_CUMVAR_CAP:-0.95}"
MEICA_SKIP_TEDANA_IF_EXISTS="${MEICA_SKIP_TEDANA_IF_EXISTS:-0}"     # 0|1

resolve_meica_priors_mat() {
  local -a candidates=()
  if [[ -n "${MEICA_PRIORS_MAT:-}" ]]; then
    candidates+=("$MEICA_PRIORS_MAT")
  fi
  if [[ -n "${NETWORK_PRIORS_MAT:-}" ]]; then
    candidates+=("$NETWORK_PRIORS_MAT")
  fi
  candidates+=("$MEDIR/res0urces/Priors.mat")
  candidates+=("$MEDIR/templates/Priors.mat")
  candidates+=("$(cd "$MEDIR/.." && pwd)/BetaPipeline/templates/Priors.mat")
  local c
  for c in "${candidates[@]}"; do
    [[ -n "$c" ]] || continue
    if [[ -f "$c" ]]; then
      readlink -f "$c"
      return 0
    fi
  done
  echo ""
}

case "$TEDANA_ACTIVATE_MODE" in
  conda_activate|conda_run|direct) ;;
  *) die "Invalid TEDANA_ACTIVATE_MODE=$TEDANA_ACTIVATE_MODE (use conda_activate|conda_run|direct)" ;;
esac
case "$TEDANA_COMPAT_MODE" in
  auto|modern|legacy|v12|v10) ;;
  *) die "Invalid TEDANA_COMPAT_MODE=$TEDANA_COMPAT_MODE (use auto|modern|legacy|v12|v10)" ;;
esac
case "$TEDANA_FITTYPE" in
  curvefit|loglin) ;;
  *) die "Invalid TEDANA_FITTYPE=$TEDANA_FITTYPE (use curvefit|loglin)" ;;
esac
case "$TEDANA_ICA_METHOD" in
  robustica|fastica) ;;
  *) die "Invalid TEDANA_ICA_METHOD=$TEDANA_ICA_METHOD (use robustica|fastica)" ;;
esac
case "$TEDANA_OVERWRITE" in
  0|1) ;;
  *) die "Invalid TEDANA_OVERWRITE=$TEDANA_OVERWRITE (use 0|1)" ;;
esac
case "$TEDANA_LOWMEM" in
  0|1) ;;
  *) die "Invalid TEDANA_LOWMEM=$TEDANA_LOWMEM (use 0|1)" ;;
esac
case "$TEDANA_USE_EXTERNAL_MIX" in
  0|1) ;;
  *) die "Invalid TEDANA_USE_EXTERNAL_MIX=$TEDANA_USE_EXTERNAL_MIX (use 0|1)" ;;
esac
case "$MEICA_RECLASSIFY_ENABLE" in
  0|1) ;;
  *) die "Invalid MEICA_RECLASSIFY_ENABLE=$MEICA_RECLASSIFY_ENABLE (use 0|1)" ;;
esac
case "$MEICA_RECLASS_NO_REPORTS" in
  0|1) ;;
  *) die "Invalid MEICA_RECLASS_NO_REPORTS=$MEICA_RECLASS_NO_REPORTS (use 0|1)" ;;
esac
case "$MEICA_CLASSIFIER_MODE" in
  nsi|legacy_template_rho|none) ;;
  tedana_rho) die "MEICA_CLASSIFIER_MODE=tedana_rho is deprecated. Use nsi or legacy_template_rho." ;;
  *) die "Invalid MEICA_CLASSIFIER_MODE=$MEICA_CLASSIFIER_MODE (use nsi|legacy_template_rho|none)" ;;
esac
case "$MEICA_QC_CIFTI_ENABLE" in
  0|1) ;;
  *) die "Invalid MEICA_QC_CIFTI_ENABLE=$MEICA_QC_CIFTI_ENABLE (use 0|1)" ;;
esac
case "$MEICA_NSI_GUARDRAIL_KAPPA_RHO" in
  0|1) ;;
  *) die "Invalid MEICA_NSI_GUARDRAIL_KAPPA_RHO=$MEICA_NSI_GUARDRAIL_KAPPA_RHO (use 0|1)" ;;
esac
case "$MEICA_SKIP_TEDANA_IF_EXISTS" in
  0|1) ;;
  *) die "Invalid MEICA_SKIP_TEDANA_IF_EXISTS=$MEICA_SKIP_TEDANA_IF_EXISTS (use 0|1)" ;;
esac
case "$MEICA_NSI_KILL_MODE" in
  adaptive|fixed) ;;
  *) die "Invalid MEICA_NSI_KILL_MODE=$MEICA_NSI_KILL_MODE (use adaptive|fixed)" ;;
esac

for cmd in python3 fslmaths fslmerge fslval awk sed parallel; do
  need_cmd "$cmd"
done
if [[ "$MEICA_QC_CIFTI_ENABLE" == "1" ]]; then
  need_cmd wb_command
fi
if [[ "$TEDANA_ACTIVATE_MODE" != "direct" ]]; then
  need_cmd conda
fi

if [[ "$MEICA_RECLASSIFY_ENABLE" == "1" && "$MEICA_CLASSIFIER_MODE" == "nsi" ]]; then
  RESOLVED_PRIORS_MAT="$(resolve_meica_priors_mat)"
  [[ -n "$RESOLVED_PRIORS_MAT" ]] || die "MEICA_CLASSIFIER_MODE=nsi but Priors.mat not found. Set NETWORK_PRIORS_MAT/MEICA_PRIORS_MAT or place priors at $MEDIR/res0urces/Priors.mat"
  if [[ "${MEICA_PRIORS_MAT:-}" != "$RESOLVED_PRIORS_MAT" ]]; then
    log "MEICA: NSI priors override -> $RESOLVED_PRIORS_MAT"
  fi
  MEICA_PRIORS_MAT="$RESOLVED_PRIORS_MAT"
fi

resolve_tedana_bin() {
  local exe="$1"
  if [[ "$TEDANA_ACTIVATE_MODE" == "direct" ]]; then
    command -v "$exe" || true
    return 0
  fi
  local base
  base="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${base:-}" && -x "$base/envs/$TEDANA_ENV/bin/$exe" ]]; then
    echo "$base/envs/$TEDANA_ENV/bin/$exe"
    return 0
  fi
  echo ""
}

TEDANA_BIN="$(resolve_tedana_bin tedana)"
if [[ -z "${TEDANA_BIN:-}" ]]; then
  die "Could not resolve tedana binary for mode=$TEDANA_ACTIVATE_MODE env=$TEDANA_ENV"
fi
ICA_RECLASSIFY_BIN="$(resolve_tedana_bin ica_reclassify)"
TEDANA_VERSION_DETECTED=""
TEDANA_CLI_PROFILE=""  # modern|legacy

tedana_exec() {
  local bin="$1"
  shift
  case "$TEDANA_ACTIVATE_MODE" in
    conda_run)
      conda run -n "$TEDANA_ENV" "$bin" "$@"
      ;;
    conda_activate)
      bash -lc '
        set -euo pipefail
        BASE="$(conda info --base 2>/dev/null || true)"
        [[ -n "${BASE:-}" && -f "$BASE/etc/profile.d/conda.sh" ]] && source "$BASE/etc/profile.d/conda.sh"
        conda activate "'"$TEDANA_ENV"'"
        exec "'"$bin"'" "$@"
      ' _ "$@"
      ;;
    direct)
      "$bin" "$@"
      ;;
  esac
}

detect_tedana_version() {
  local raw ver
  raw="$(tedana_exec "$TEDANA_BIN" --version 2>&1 | head -n 1 || true)"
  if [[ -z "${raw// }" ]]; then
    raw="$(tedana_exec "$TEDANA_BIN" -h 2>&1 | head -n 5 || true)"
  fi
  ver="$(echo "$raw" | grep -Eo '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n 1 || true)"
  echo "$ver"
}

resolve_tedana_cli_profile() {
  local mode="$TEDANA_COMPAT_MODE"
  case "$mode" in
    modern|v12)
      TEDANA_CLI_PROFILE="modern"
      TEDANA_VERSION_DETECTED="$(detect_tedana_version)"
      ;;
    legacy|v10)
      TEDANA_CLI_PROFILE="legacy"
      TEDANA_VERSION_DETECTED="$(detect_tedana_version)"
      ;;
    auto)
      TEDANA_VERSION_DETECTED="$(detect_tedana_version)"
      if [[ -z "$TEDANA_VERSION_DETECTED" ]]; then
        log "MEICA: could not parse tedana version, defaulting CLI profile to modern"
        TEDANA_CLI_PROFILE="modern"
        return 0
      fi
      local major="${TEDANA_VERSION_DETECTED%%.*}"
      if [[ "$major" =~ ^[0-9]+$ ]] && [[ "$major" -eq 0 ]]; then
        TEDANA_CLI_PROFILE="legacy"
      else
        TEDANA_CLI_PROFILE="modern"
      fi
      ;;
    *)
      die "Unhandled TEDANA_COMPAT_MODE=$mode"
      ;;
  esac
}

sanity_check_tedana() {
  resolve_tedana_cli_profile
  log "MEICA: tedana mode=$TEDANA_ACTIVATE_MODE env=$TEDANA_ENV compat_request=$TEDANA_COMPAT_MODE detected_version=${TEDANA_VERSION_DETECTED:-unknown} cli_profile=$TEDANA_CLI_PROFILE bin=$TEDANA_BIN"
  if ! tedana_exec "$TEDANA_BIN" -h >/dev/null 2>&1; then
    die "Tedana sanity check failed for env '$TEDANA_ENV'. Install a working tedana in that env (or use TEDANA_ACTIVATE_MODE=direct), then retry."
  fi
  if [[ "$MEICA_RECLASSIFY_ENABLE" == "1" && "$TEDANA_CLI_PROFILE" != "legacy" ]]; then
    [[ -x "${ICA_RECLASSIFY_BIN:-}" ]] || die "Could not resolve ica_reclassify binary for env '$TEDANA_ENV'"
    [[ -f "$MEICA_RECLASSIFY_PY" ]] || die "Missing classifier script: $MEICA_RECLASSIFY_PY"
    if ! tedana_exec "$ICA_RECLASSIFY_BIN" -h >/dev/null 2>&1; then
      die "ica_reclassify sanity check failed for env '$TEDANA_ENV'"
    fi
  fi
}

list_runs() {
  local subdir="$1"
  local start="$2"
  mapfile -t SESS < <(find "$subdir/func/$FuncDirName" -maxdepth 1 -type d -name 'session_*' | sort)
  for sesdir in "${SESS[@]}"; do
    local s="${sesdir##*/}"; s="${s#session_}"
    [[ "$s" -ge "$start" ]] || continue
    find "$sesdir" -maxdepth 1 -type d -name 'run_*' | sort
  done
}

discover_echoes_sorted() {
  local rundir="$1"
  python3 - "$rundir" "$FuncFilePrefix" <<'PY'
import glob, os, re, sys
rundir = sys.argv[1]
prefix = sys.argv[2]
rx = re.compile(re.escape(prefix) + r"_E(\d+)_acpc\.nii\.gz$")
pairs = []
for f in glob.glob(os.path.join(rundir, f"{prefix}_E*_acpc.nii.gz")):
    m = rx.search(os.path.basename(f))
    if m:
        pairs.append((int(m.group(1)), f))
pairs.sort()
for _, f in pairs:
    print(f)
PY
}

te_args_from_file() {
  local tefile="$1"
  local cli_profile="${2:-modern}" # modern|legacy
  [[ -f "$tefile" ]] || return 1
  python3 - "$tefile" "$cli_profile" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
cli_profile = str(sys.argv[2]).strip().lower()
vals = []
for tok in p.read_text().split():
    vals.append(float(tok))
if not vals:
    raise SystemExit("No TE values found")

# Legacy tedana v0.x behavior expects TE in milliseconds.
# Modern tedana (>=1.x) expects TE in seconds.
if cli_profile == "legacy":
    # If file appears to be in seconds, convert to ms for legacy mode.
    if all(v <= 1.0 for v in vals):
        vals = [v * 1000.0 for v in vals]
else:
    # If file appears to be in ms, convert to seconds for modern mode.
    if any(v > 1.0 for v in vals):
        vals = [v / 1000.0 for v in vals]

if any(vals[i] <= vals[i - 1] for i in range(1, len(vals))):
    raise SystemExit("TE values must be strictly ascending")

for v in vals:
    print(f"{v:.8f}".rstrip("0").rstrip("."))
PY
}

map_qc_volume_to_cifti() {
  local in_vol="$1"
  local out_cifti="$2"
  local tmpdir="$3"
  local dim4="$4"

  local pial_l="$Subdir/anat/T1w/Native/${Subject}.L.pial.native.surf.gii"
  local white_l="$Subdir/anat/T1w/Native/${Subject}.L.white.native.surf.gii"
  local mid_l="$Subdir/anat/T1w/Native/${Subject}.L.midthickness.native.surf.gii"
  local pial_r="$Subdir/anat/T1w/Native/${Subject}.R.pial.native.surf.gii"
  local white_r="$Subdir/anat/T1w/Native/${Subject}.R.white.native.surf.gii"
  local mid_r="$Subdir/anat/T1w/Native/${Subject}.R.midthickness.native.surf.gii"

  local roi_l="$Subdir/anat/MNINonLinear/Native/${Subject}.L.roi.native.shape.gii"
  local roi_r="$Subdir/anat/MNINonLinear/Native/${Subject}.R.roi.native.shape.gii"
  local roi32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.atlasroi.32k_fs_LR.shape.gii"
  local roi32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.atlasroi.32k_fs_LR.shape.gii"
  local reg_l="$Subdir/anat/MNINonLinear/Native/${Subject}.L.sphere.MSMSulc.native.surf.gii"
  local reg_r="$Subdir/anat/MNINonLinear/Native/${Subject}.R.sphere.MSMSulc.native.surf.gii"
  local reg32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.sphere.32k_fs_LR.surf.gii"
  local reg32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.sphere.32k_fs_LR.surf.gii"
  local mid32_l="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.L.midthickness.32k_fs_LR.surf.gii"
  local mid32_r="$Subdir/anat/MNINonLinear/fsaverage_LR32k/${Subject}.R.midthickness.32k_fs_LR.surf.gii"
  local subcort="$Subdir/func/rois/Subcortical_ROIs_acpc.nii.gz"

  for f in "$pial_l" "$white_l" "$mid_l" "$pial_r" "$white_r" "$mid_r" \
           "$roi_l" "$roi_r" "$roi32_l" "$roi32_r" "$reg_l" "$reg_r" \
           "$reg32_l" "$reg32_r" "$mid32_l" "$mid32_r" "$subcort"; do
    [[ -f "$f" ]] || die "MEICA QC CIFTI: missing required file: $f"
  done

  rm -rf "$tmpdir" 2>/dev/null || true
  mkdir -p "$tmpdir"
  local l_native="$tmpdir/lh.native.shape.gii"
  local r_native="$tmpdir/rh.native.shape.gii"
  local l_32k="$tmpdir/lh.32k_fs_LR.shape.gii"
  local r_32k="$tmpdir/rh.32k_fs_LR.shape.gii"

  wb_command -volume-to-surface-mapping "$in_vol" "$mid_l" "$l_native" -ribbon-constrained "$white_l" "$pial_l"
  wb_command -volume-to-surface-mapping "$in_vol" "$mid_r" "$r_native" -ribbon-constrained "$white_r" "$pial_r"
  wb_command -metric-dilate "$l_native" "$mid_l" 10 "$l_native" -nearest
  wb_command -metric-dilate "$r_native" "$mid_r" 10 "$r_native" -nearest
  wb_command -metric-mask "$l_native" "$roi_l" "$l_native"
  wb_command -metric-mask "$r_native" "$roi_r" "$r_native"
  wb_command -metric-resample "$l_native" "$reg_l" "$reg32_l" ADAP_BARY_AREA "$l_32k" -area-surfs "$mid_l" "$mid32_l" -current-roi "$roi_l"
  wb_command -metric-resample "$r_native" "$reg_r" "$reg32_r" ADAP_BARY_AREA "$r_32k" -area-surfs "$mid_r" "$mid32_r" -current-roi "$roi_r"
  wb_command -metric-mask "$l_32k" "$roi32_l" "$l_32k"
  wb_command -metric-mask "$r_32k" "$roi32_r" "$r_32k"

  if [[ "$dim4" -gt 1 ]]; then
    wb_command -cifti-create-dense-timeseries "$out_cifti" \
      -volume "$in_vol" "$subcort" \
      -left-metric "$l_32k" -roi-left "$roi32_l" \
      -right-metric "$r_32k" -roi-right "$roi32_r"
  else
    wb_command -cifti-create-dense-scalar "$out_cifti" \
      -volume "$in_vol" "$subcort" \
      -left-metric "$l_32k" \
      -right-metric "$r_32k"
  fi
  rm -rf "$tmpdir" 2>/dev/null || true
}

create_orig_aliases_from_modern() {
  local teddir="$1"
  [[ "$MEICA_ORIG_ALIAS_ENABLE" == "1" ]] || return 0

  # tedana >=1.x may emit desc-* outputs regardless of --convention; expose
  # legacy/orig aliases expected by downstream reclassify/QC steps.
  local -a pairs=(
    "desc-ICA_stat-z_components.nii.gz:betas_OC.nii.gz"
    "desc-ICAAccepted_components.nii.gz:betas_hik_OC.nii.gz"
    "desc-ICAAccepted_stat-z_components.nii.gz:feats_OC2.nii.gz"
    "desc-denoised_bold.nii.gz:dn_ts_OC.nii.gz"
    "desc-optcom_bold.nii.gz:ts_OC.nii.gz"
    "desc-ICA_mixing.tsv:ica_mixing.tsv"
    "desc-ICA_decomposition.json:ica_decomposition.json"
    "desc-tedana_metrics.tsv:ica_metrics.tsv"
    "desc-tedana_metrics.json:ica_metrics.json"
    "desc-ICACrossComponent_metrics.json:ica_cross_component_metrics.json"
    "desc-ICA_status_table.tsv:ica_status_table.tsv"
    "T2starmap.nii.gz:t2svG.nii.gz"
    "S0map.nii.gz:s0vG.nii.gz"
  )

  local pair src dst
  for pair in "${pairs[@]}"; do
    src="${pair%%:*}"
    dst="${pair##*:}"
    if [[ -f "$teddir/$src" && ! -e "$teddir/$dst" ]]; then
      ln -s "$src" "$teddir/$dst"
      log "MEICA: created orig alias $dst -> $src"
    fi
  done
}

generate_meica_qc_ciftis() {
  local rundir="$1"
  local teddir="$2"
  local tags_csv="$3"

  IFS=',' read -r -a tags <<<"$tags_csv"
  for raw_tag in "${tags[@]}"; do
    local tag
    tag="$(echo "$raw_tag" | sed 's/^ *//; s/ *$//')"
    [[ -n "$tag" ]] || continue

    local in_vol=""
    local out_cifti=""
    case "$tag" in
      betas_OC)
        in_vol="$teddir/betas_OC.nii.gz"
        out_cifti="$teddir/betas_OC.dtseries.nii"
        ;;
      t2sv)
        if [[ -f "$teddir/t2sv.nii.gz" ]]; then
          in_vol="$teddir/t2sv.nii.gz"
        else
          in_vol="$teddir/t2svG.nii.gz"
        fi
        out_cifti="$teddir/t2sv.dscalar.nii"
        ;;
      s0v)
        if [[ -f "$teddir/s0v.nii.gz" ]]; then
          in_vol="$teddir/s0v.nii.gz"
        else
          in_vol="$teddir/s0vG.nii.gz"
        fi
        out_cifti="$teddir/s0v.dscalar.nii"
        ;;
      *)
        log "MEICA QC CIFTI: skipping unknown tag '$tag'"
        continue
        ;;
    esac

    if [[ ! -f "$in_vol" ]]; then
      log "MEICA QC CIFTI: input missing for $tag: $in_vol (skip)"
      continue
    fi
    local dim4
    dim4="$(fslval "$in_vol" dim4 2>/dev/null || echo 1)"
    [[ -n "$dim4" ]] || dim4=1
    map_qc_volume_to_cifti "$in_vol" "$out_cifti" "$rundir/.tmp_meica_cifti_${tag}" "$dim4"
    log "MEICA QC CIFTI: wrote $(basename "$out_cifti")"
  done
}

make_component_subset_nifti() {
  local in_4d="$1"
  local ids_file="$2"
  local out_4d="$3"
  python3 - "$in_4d" "$ids_file" "$out_4d" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np

in_4d = Path(sys.argv[1])
ids_file = Path(sys.argv[2])
out_4d = Path(sys.argv[3])

img = nib.load(str(in_4d))
data = np.asarray(img.dataobj, dtype=np.float64)
if data.ndim != 4:
    raise SystemExit(f"Expected 4D NIfTI: {in_4d}")

ids = []
if ids_file.exists():
    for tok in ids_file.read_text().split():
        ids.append(int(tok))
ids = [i for i in ids if 0 <= i < data.shape[3]]

if len(ids) == 0:
    subset = np.zeros(data.shape[:3] + (1,), dtype=np.float64)
else:
    subset = data[:, :, :, ids]

nib.save(nib.Nifti1Image(subset, img.affine, img.header), str(out_4d))
print(len(ids))
PY
}

extract_component_ids_from_decisions() {
  local decisions_tsv="$1"
  local column_name="$2"
  local out_ids="$3"
  python3 - "$decisions_tsv" "$column_name" "$out_ids" <<'PY'
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t")
col = sys.argv[2]
out = sys.argv[3]
if col not in df.columns:
    open(out, "w").close()
    raise SystemExit(0)
ids = df.loc[df[col].astype(bool), "component_id"].astype(int).tolist()
with open(out, "w") as f:
    f.write("\n".join(str(i) for i in ids))
    if ids:
        f.write("\n")
PY
}

make_reclass_qc_ciftis() {
  local rundir="$1"
  local reclass_dir="$2"
  local src_tedana="$3"

  local base_vol="$reclass_dir/betas_OC.nii.gz"
  [[ -f "$base_vol" ]] || return 0

  local full_dt="$reclass_dir/betas_OC.dtseries.nii"
  map_qc_volume_to_cifti "$base_vol" "$full_dt" "$rundir/.tmp_meica_cifti_reclass_full" "$(fslval "$base_vol" dim4)"

  local acc_ids="$reclass_dir/AcceptedComponents.txt"
  local rej_ids="$reclass_dir/RejectedComponents.txt"
  local acc_vol="$reclass_dir/betas_OC_Accepted.nii.gz"
  local rej_vol="$reclass_dir/betas_OC_Rejected.nii.gz"
  make_component_subset_nifti "$base_vol" "$acc_ids" "$acc_vol" >/dev/null
  make_component_subset_nifti "$base_vol" "$rej_ids" "$rej_vol" >/dev/null
  map_qc_volume_to_cifti "$acc_vol" "$reclass_dir/betas_OC_Accepted.dtseries.nii" "$rundir/.tmp_meica_cifti_reclass_acc" "$(fslval "$acc_vol" dim4)"
  map_qc_volume_to_cifti "$rej_vol" "$reclass_dir/betas_OC_Rejected.dtseries.nii" "$rundir/.tmp_meica_cifti_reclass_rej" "$(fslval "$rej_vol" dim4)"

  local decisions="$reclass_dir/ComponentDecisions.tsv"
  if [[ -f "$decisions" ]]; then
    local rescued_ids="$rundir/.tmp_meica_rescued_ids.txt"
    local killed_ids="$rundir/.tmp_meica_killed_ids.txt"
    extract_component_ids_from_decisions "$decisions" "rescued" "$rescued_ids"
    extract_component_ids_from_decisions "$decisions" "killed" "$killed_ids"
    local rescued_vol="$reclass_dir/betas_OC_Rescued.nii.gz"
    local killed_vol="$reclass_dir/betas_OC_Killed.nii.gz"
    make_component_subset_nifti "$base_vol" "$rescued_ids" "$rescued_vol" >/dev/null
    make_component_subset_nifti "$base_vol" "$killed_ids" "$killed_vol" >/dev/null
    map_qc_volume_to_cifti "$rescued_vol" "$reclass_dir/betas_OC_Rescued.dtseries.nii" "$rundir/.tmp_meica_cifti_reclass_resc" "$(fslval "$rescued_vol" dim4)"
    map_qc_volume_to_cifti "$killed_vol" "$reclass_dir/betas_OC_Killed.dtseries.nii" "$rundir/.tmp_meica_cifti_reclass_kill" "$(fslval "$killed_vol" dim4)"
    rm -f "$rescued_ids" "$killed_ids"
  fi
}

process_run() {
  local rundir="$1"
  local rel="${rundir#$Subdir/func/$FuncDirName/}"
  local outdir="$rundir/$MEICA_TEDANA_SUBDIR"
  local mask="$rundir/brain_mask.nii.gz"
  local tmp="$rundir/tmp.nii.gz"
  local reclass_dir="$outdir/$MEICA_RECLASS_SUBDIR"
  cleanup_run_tmp() {
    rm -f "$mask" "$tmp"
  }
  trap cleanup_run_tmp RETURN
  local has_tedana_outputs=0
  if [[ "$TEDANA_CLI_PROFILE" == "legacy" ]]; then
    if [[ -f "$outdir/ica_decomposition.json" && -f "$outdir/ica_mixing.tsv" && -f "$outdir/dn_ts_OC.nii.gz" ]]; then
      has_tedana_outputs=1
    fi
  else
    if [[ ( -f "$outdir/registry.json" || -f "$outdir/desc-tedana_registry.json" ) && \
          ( -f "$outdir/ica_metrics.tsv" || -f "$outdir/desc-tedana_metrics.tsv" ) ]]; then
      has_tedana_outputs=1
    fi
  fi

  log "MEICA: start $rel"
  if [[ "$MEICA_SKIP_TEDANA_IF_EXISTS" == "1" && "$has_tedana_outputs" == "1" ]]; then
    log "MEICA: preserving existing Tedana folder for $rel (MEICA_SKIP_TEDANA_IF_EXISTS=1)"
  else
    rm -rf "$outdir" 2>/dev/null || true
    mkdir -p "$outdir"
  fi
  rm -rf "$reclass_dir" "$rundir/Tedana+ManualComponentClassification" 2>/dev/null || true
  mkdir -p "$outdir"

  # Build a conservative binary mask from anatomical brain in functional space and signal-positive voxels.
  fslmaths "$rundir/${FuncFilePrefix}_E1_acpc.nii.gz" -Tmin "$tmp"
  fslmaths "$Subdir/func/xfms/$FuncXfmsDir/T1w_acpc_brain_func.nii.gz" -mas "$tmp" -bin "$mask"

  mapfile -t ECHOES < <(discover_echoes_sorted "$rundir")
  [[ "${#ECHOES[@]}" -gt 0 ]] || die "MEICA: no ${FuncFilePrefix}_E*_acpc.nii.gz files in $rundir"

  mapfile -t TE_ARR < <(te_args_from_file "$rundir/TE.txt" "$TEDANA_CLI_PROFILE")
  [[ "${#TE_ARR[@]}" -gt 0 ]] || die "MEICA: failed to parse TE list in $rundir/TE.txt"

  if [[ "${#ECHOES[@]}" -ne "${#TE_ARR[@]}" ]]; then
    die "MEICA: TE count (${#TE_ARR[@]}) != echo count (${#ECHOES[@]}) in $rundir"
  fi

  log "MEICA: tedana run $rel (echoes=${#ECHOES[@]} env=$TEDANA_ENV mode=$TEDANA_ACTIVATE_MODE)"
  local -a TEDANA_ARGS=()
  TEDANA_ARGS+=( -d "${ECHOES[@]}" )
  TEDANA_ARGS+=( -e "${TE_ARR[@]}" )
  TEDANA_ARGS+=( --out-dir "$outdir" )
  TEDANA_ARGS+=( --mask "$mask" )
  if [[ "$TEDANA_CLI_PROFILE" == "modern" ]]; then
    TEDANA_ARGS+=( --prefix "" )
    TEDANA_ARGS+=( --convention "$TEDANA_CONVENTION" )
    TEDANA_ARGS+=( --masktype "$TEDANA_MASKTYPE" )
  fi
  TEDANA_ARGS+=( --fittype "$TEDANA_FITTYPE" )
  TEDANA_ARGS+=( --combmode t2s )
  TEDANA_ARGS+=( --tedpca "$MEPCA" )
  if [[ "$TEDANA_CLI_PROFILE" == "modern" ]]; then
    TEDANA_ARGS+=( --ica-method "$TEDANA_ICA_METHOD" )
  fi
  TEDANA_ARGS+=( --seed "$TEDANA_SEED" )
  TEDANA_ARGS+=( --maxit "$MaxIterations" )
  TEDANA_ARGS+=( --maxrestart "$MaxRestarts" )
  TEDANA_ARGS+=( --n-threads "$TEDANA_THREADS" )
  if [[ "$TEDANA_CLI_PROFILE" == "modern" && "$TEDANA_ICA_METHOD" == "robustica" ]]; then
    TEDANA_ARGS+=( --n-robust-runs "$TEDANA_N_ROBUST_RUNS" )
  fi
  if [[ "$TEDANA_LOWMEM" == "1" ]]; then
    TEDANA_ARGS+=( --lowmem )
  fi
  if [[ "$TEDANA_CLI_PROFILE" == "modern" && "$TEDANA_OVERWRITE" == "1" ]]; then
    TEDANA_ARGS+=( --overwrite )
  fi
  if [[ "$TEDANA_CLI_PROFILE" == "modern" && "$TEDANA_USE_EXTERNAL_MIX" == "1" ]]; then
    [[ -n "$TEDANA_EXTERNAL_MIX_BASENAME" ]] || die "TEDANA_USE_EXTERNAL_MIX=1 requires TEDANA_EXTERNAL_MIX_BASENAME"
    local ext_mix="$rundir/$TEDANA_EXTERNAL_MIX_BASENAME"
    [[ -f "$ext_mix" ]] || die "Missing external tedana mix for $rel: $ext_mix"
    TEDANA_ARGS+=( --mix "$ext_mix" )
    log "MEICA: using external mix for $rel -> $ext_mix"
  fi

  if [[ "$MEICA_SKIP_TEDANA_IF_EXISTS" == "1" && "$has_tedana_outputs" == "1" ]]; then
    log "MEICA: skipping tedana execution for $rel (existing outputs detected)"
  else
    tedana_exec "$TEDANA_BIN" "${TEDANA_ARGS[@]}"
  fi

  create_orig_aliases_from_modern "$outdir"

  if [[ "$MEICA_QC_CIFTI_ENABLE" == "1" ]]; then
    generate_meica_qc_ciftis "$rundir" "$outdir" "$MEICA_QC_CIFTI_TAGS"
  fi

  local final_dir="$outdir"
  if [[ "$MEICA_RECLASSIFY_ENABLE" == "1" ]]; then
    log "MEICA: auto-screen + reclassify $rel (mode=$MEICA_CLASSIFIER_MODE)"
    local -a AUTO_ARGS=()
    AUTO_ARGS+=( --data-dir "$rundir" )
    AUTO_ARGS+=( --tedana-dir "$outdir" )
    AUTO_ARGS+=( --out-dir "$reclass_dir" )
    AUTO_ARGS+=( --classifier-mode "$MEICA_CLASSIFIER_MODE" )
    AUTO_ARGS+=( --rho-rescue "$MEICA_RHO_RESCUE" )
    AUTO_ARGS+=( --rho-reject "$MEICA_RHO_REJECT" )
    AUTO_ARGS+=( --rescue-nsi "$MEICA_NSI_RESCUE_THRESHOLD" )
    AUTO_ARGS+=( --rescue-quantile "$MEICA_NSI_RESCUE_QUANTILE" )
    AUTO_ARGS+=( --kill-mode "$MEICA_NSI_KILL_MODE" )
    AUTO_ARGS+=( --kill-nsi "$MEICA_NSI_KILL_THRESHOLD" )
    AUTO_ARGS+=( --kill-nsi-min "$MEICA_NSI_KILL_MIN" )
    AUTO_ARGS+=( --kill-nsi-max "$MEICA_NSI_KILL_MAX" )
    AUTO_ARGS+=( --kill-intercept "$MEICA_NSI_KILL_INTERCEPT" )
    AUTO_ARGS+=( --kill-slope "$MEICA_NSI_KILL_SLOPE" )
    AUTO_ARGS+=( --guardrail-kappa-rho "$MEICA_NSI_GUARDRAIL_KAPPA_RHO" )
    AUTO_ARGS+=( --subcort-ratio-thresh "$MEICA_SUBCORT_RATIO_THRESH" )
    AUTO_ARGS+=( --kill-priority-enable "$MEICA_KILL_PRIORITY_ENABLE" )
    AUTO_ARGS+=( --kill-priority-w-logratio "$MEICA_KILL_PRIORITY_W_LOGRATIO" )
    AUTO_ARGS+=( --kill-priority-w-nsi "$MEICA_KILL_PRIORITY_W_NSI" )
    AUTO_ARGS+=( --kill-priority-w-var "$MEICA_KILL_PRIORITY_W_VAR" )
    AUTO_ARGS+=( --kill-var-floor-quantile "$MEICA_KILL_VAR_FLOOR_QUANTILE" )
    AUTO_ARGS+=( --kill-cumvar-cap "$MEICA_KILL_CUMVAR_CAP" )
    if [[ -n "${MEICA_PRIORS_MAT:-}" ]]; then
      AUTO_ARGS+=( --priors-mat "$MEICA_PRIORS_MAT" )
    fi
    local betas_cifti="${MEICA_BETAS_CIFTI:-}"
    if [[ -z "$betas_cifti" && -f "$outdir/betas_OC.dtseries.nii" ]]; then
      betas_cifti="$outdir/betas_OC.dtseries.nii"
    fi
    if [[ -n "$betas_cifti" ]]; then
      AUTO_ARGS+=( --betas-cifti "$betas_cifti" )
    fi
    python3 "$MEICA_RECLASSIFY_PY" "${AUTO_ARGS[@]}"

    local acc="$reclass_dir/AcceptedComponents.txt"
    local rej="$reclass_dir/RejectedComponents.txt"
    [[ -f "$acc" ]] || die "MEICA: classifier missing $acc"
    [[ -f "$rej" ]] || die "MEICA: classifier missing $rej"

    if [[ "$TEDANA_CLI_PROFILE" == "legacy" ]]; then
      local mix=""
      local ctab=""
      for p in "$outdir/ica_mixing.tsv" "$outdir/ica_mix.tsv"; do
        if [[ -f "$p" ]]; then mix="$p"; break; fi
      done
      for p in "$outdir/ica_decomposition.json" "$outdir/ica_metrics.json"; do
        if [[ -f "$p" ]]; then ctab="$p"; break; fi
      done
      [[ -n "$mix" ]] || die "MEICA v10 mode: missing ICA mixing TSV in $outdir"
      [[ -n "$ctab" ]] || die "MEICA v10 mode: missing decomposition JSON in $outdir"

      mapfile -t MANACC_IDS < <(awk 'NF{print $1}' "$acc")
      [[ "${#MANACC_IDS[@]}" -gt 0 ]] || die "MEICA v10 mode: AcceptedComponents.txt is empty; cannot pass --manacc"

      mkdir -p "$reclass_dir"
      local -a RECLASS_TEDANA_ARGS=()
      RECLASS_TEDANA_ARGS+=( -d "${ECHOES[@]}" )
      RECLASS_TEDANA_ARGS+=( -e "${TE_ARR[@]}" )
      RECLASS_TEDANA_ARGS+=( --out-dir "$reclass_dir" )
      RECLASS_TEDANA_ARGS+=( --mask "$mask" )
      RECLASS_TEDANA_ARGS+=( --fittype "$TEDANA_FITTYPE" )
      RECLASS_TEDANA_ARGS+=( --combmode t2s )
      RECLASS_TEDANA_ARGS+=( --tedpca "$MEPCA" )
      RECLASS_TEDANA_ARGS+=( --seed "$TEDANA_SEED" )
      RECLASS_TEDANA_ARGS+=( --maxit "$MaxIterations" )
      RECLASS_TEDANA_ARGS+=( --maxrestart "$MaxRestarts" )
      RECLASS_TEDANA_ARGS+=( --n-threads "$TEDANA_THREADS" )
      RECLASS_TEDANA_ARGS+=( --mix "$mix" )
      RECLASS_TEDANA_ARGS+=( --ctab "$ctab" )
      RECLASS_TEDANA_ARGS+=( --manacc "${MANACC_IDS[@]}" )
      if [[ "$TEDANA_LOWMEM" == "1" ]]; then
        RECLASS_TEDANA_ARGS+=( --lowmem )
      fi
      if [[ "$MEICA_RECLASS_NO_REPORTS" == "1" ]]; then
        RECLASS_TEDANA_ARGS+=( --no-reports )
      fi
      tedana_exec "$TEDANA_BIN" "${RECLASS_TEDANA_ARGS[@]}"
    else
      local registry="$outdir/registry.json"
      if [[ ! -f "$registry" ]]; then
        registry="$outdir/desc-tedana_registry.json"
      fi
      [[ -f "$registry" ]] || die "MEICA: missing registry.json/desc-tedana_registry.json for reclassify in $outdir"
      local -a RECLASS_ARGS=()
      RECLASS_ARGS+=( "$registry" )
      RECLASS_ARGS+=( --manacc "$acc" )
      RECLASS_ARGS+=( --manrej "$rej" )
      RECLASS_ARGS+=( --out-dir "$reclass_dir" )
      RECLASS_ARGS+=( --convention "$TEDANA_CONVENTION" )
      RECLASS_ARGS+=( --prefix "" )
      if [[ "$MEICA_RECLASS_NO_REPORTS" == "1" ]]; then
        RECLASS_ARGS+=( --no-reports )
      fi
      if [[ "$TEDANA_OVERWRITE" == "1" ]]; then
        RECLASS_ARGS+=( --overwrite )
      fi
      tedana_exec "$ICA_RECLASSIFY_BIN" "${RECLASS_ARGS[@]}"
    fi
    if [[ "$MEICA_QC_CIFTI_ENABLE" == "1" ]]; then
      make_reclass_qc_ciftis "$rundir" "$reclass_dir" "$outdir"
      log "MEICA QC CIFTI: wrote reclass component CIFTIs in $reclass_dir"
    fi
    final_dir="$reclass_dir"
  fi

  # Map tedana v11 outputs to legacy pipeline filenames for downstream compatibility.
  local oc=""
  local dn=""
  for c in \
    "$final_dir/desc-optcom_bold.nii.gz" \
    "$final_dir/ts_OC.nii.gz" \
    "$outdir/desc-optcom_bold.nii.gz" \
    "$outdir/ts_OC.nii.gz"; do
    if [[ -f "$c" ]]; then oc="$c"; break; fi
  done
  for c in \
    "$final_dir/desc-denoised_bold.nii.gz" \
    "$final_dir/dn_ts_OC.nii.gz" \
    "$outdir/desc-denoised_bold.nii.gz" \
    "$outdir/dn_ts_OC.nii.gz"; do
    if [[ -f "$c" ]]; then dn="$c"; break; fi
  done
  [[ -n "$oc" ]] || die "MEICA: tedana/reclass output missing optcom time-series in $final_dir"
  [[ -n "$dn" ]] || die "MEICA: tedana/reclass output missing denoised time-series in $final_dir"
  cp -f "$oc" "$rundir/${FuncFilePrefix}_OCME.nii.gz"
  cp -f "$dn" "$rundir/${FuncFilePrefix}_OCME+MEICA.nii.gz"

  log "MEICA: done $rel"
}

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
if [[ "$WORKER_MODE" == "1" ]]; then
  [[ -n "$WORKER_RUNDIR" ]] || die "MEICA worker mode missing run directory argument"
  rundir="$WORKER_RUNDIR"
  Subdir="$StudyFolder/$Subject"
  process_run "$rundir"
  exit 0
fi

log "MEICA: start subject=$Subject StartSession=$StartSession NTHREADS=$NTHREADS jobs=$MEICA_PARALLEL_JOBS"
log "MEICA: tedana env=$TEDANA_ENV mode=$TEDANA_ACTIVATE_MODE compat_request=$TEDANA_COMPAT_MODE fittype=$TEDANA_FITTYPE ica=$TEDANA_ICA_METHOD"
log "MEICA: functional naming func/$FuncDirName with prefix ${FuncFilePrefix}_*"

mapfile -t RUNS < <(list_runs "$Subdir" "$StartSession")
[[ "${#RUNS[@]}" -gt 0 ]] || die "MEICA: no run directories found in $Subdir/func/$FuncDirName"
sanity_check_tedana

if [[ "$MEICA_PARALLEL_JOBS" -le 1 ]]; then
  for r in "${RUNS[@]}"; do
    process_run "$r"
  done
else
  printf "%s\n" "${RUNS[@]}" | parallel --jobs "$MEICA_PARALLEL_JOBS" \
    "$SELF" --worker "$Subject" "$StudyFolder" "$NTHREADS" "$MEPCA" "$MaxIterations" "$MaxRestarts" "$StartSession" "$MEDIR" {}
fi

log "MEICA: all done subject=$Subject"

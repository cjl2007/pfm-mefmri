#!/usr/bin/env bash
# run_charm_and_masks.sh
# CJL; (cjl2007@med.cornell.edu)
#
# Purpose:
#   1) Run SimNIBS CHARM using HCP outputs in ${StudyFolder}/${Subject}/anat/T1w
#   2) Rename m2m_${Subject} -> anat/charm and stash log
#   3) Cleanup aparc+aseg using CHARM labeling with targeted ROI retention
#   4) Create a binarized brain mask from CHARM labeling (labeling < 100)
#   5) Apply that mask to:
#        - T1w_acpc.nii.gz               -> T1w_acpc_brain.nii.gz
#        - T1w_acpc_dc_restore.nii.gz    -> T1w_acpc_dc_restore_brain.nii.gz
#        - T2w_acpc_dc_restore.nii.gz    -> T2w_acpc_dc_restore_brain.nii.gz   (ADDED)
#

set -euo pipefail

StudyFolder="${1:?Need StudyFolder}"
Subject="${2:?Need Subject}"
CHARM_BIN="${3:-${CHARM_BIN:-charm}}"

# ---- paths ----
ANAT_DIR="${StudyFolder}/${Subject}/anat"
T1W_DIR="${ANAT_DIR}/T1w"

CHARM_T1="${T1W_DIR}/T1w_acpc_dc_restore.nii.gz"
CHARM_T2="${T1W_DIR}/T2w_acpc_dc_restore.nii.gz"

# CHARM outputs
CHARM_TISSUE="${ANAT_DIR}/charm/final_tissues.nii.gz"
CHARM_LABELING="${ANAT_DIR}/charm/segmentation/labeling.nii.gz"

# HCP/Freesurfer aparc+aseg moved into anat/T1w
ASEG_FILE="${T1W_DIR}/aparc+aseg.nii.gz"

# Output brain mask (in same directory as brainMask.nii.gz would live)
BRAIN_MASK_OUT="${T1W_DIR}/T1w_acpc_brain_mask.nii.gz"

# Images to mask + outputs
T1W_ACPC="${T1W_DIR}/T1w_acpc.nii.gz"
T1W_ACPC_BRAIN="${T1W_DIR}/T1w_acpc_brain.nii.gz"

T1W_ACPC_DCRESTORE="${T1W_DIR}/T1w_acpc_dc_restore.nii.gz"
T1W_ACPC_DCRESTORE_BRAIN="${T1W_DIR}/T1w_acpc_dc_restore_brain.nii.gz"

T2W_ACPC_DCRESTORE="${T1W_DIR}/T2w_acpc_dc_restore.nii.gz"
T2W_ACPC_DCRESTORE_BRAIN="${T1W_DIR}/T2w_acpc_dc_restore_brain.nii.gz"

# Minimal logging to file
LOG_FILE="${ANAT_DIR}/charm.txt"
CHARM_REUSE_EXISTING_M2M="${CHARM_REUSE_EXISTING_M2M:-0}" # 0|1
CHARM_SKIP_RUN="${CHARM_SKIP_RUN:-0}" # 0|1; 1 => skip CHARM execution and use existing anat/charm outputs
CHARM_BRAIN_MASK_MODE="${CHARM_BRAIN_MASK_MODE:-charm}" # charm|hcp
CHARM_BRAIN_MASK_DILATE_ITERS="${CHARM_BRAIN_MASK_DILATE_ITERS:-1}" # integer >= 0 (used in charm mode)
CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS="${CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS:-1}" # 0|1
CHARM_WRITE_CORTICAL_RIBBON="${CHARM_WRITE_CORTICAL_RIBBON:-1}" # 0|1

# ---- sanity checks ----
if [ ! -d "${ANAT_DIR}" ]; then
  echo "ERROR: Missing ANAT_DIR: ${ANAT_DIR}" >&2
  exit 1
fi

if [ ! -d "${T1W_DIR}" ]; then
  echo "ERROR: Missing T1W_DIR: ${T1W_DIR}" >&2
  exit 1
fi

if [ "${CHARM_SKIP_RUN}" != "1" ]; then
  if [[ "${CHARM_BIN}" == */* ]]; then
    if [ ! -x "${CHARM_BIN}" ]; then
      echo "ERROR: CHARM binary not found or not executable: ${CHARM_BIN}" >&2
      exit 1
    fi
  else
    if ! CHARM_BIN_RESOLVED="$(command -v "${CHARM_BIN}" 2>/dev/null)"; then
      echo "ERROR: CHARM command not found on PATH: ${CHARM_BIN}" >&2
      echo "Set CHARM_BIN in your config to an absolute path, or add CHARM to PATH." >&2
      exit 1
    fi
    CHARM_BIN="${CHARM_BIN_RESOLVED}"
  fi
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 not found in PATH (needed for CHARM socket preflight)." >&2
  exit 1
fi

if [ ! -f "${CHARM_T1}" ]; then
  echo "ERROR: Missing CHARM T1 input: ${CHARM_T1}" >&2
  exit 1
fi

if ! command -v fslmaths >/dev/null 2>&1; then
  echo "ERROR: fslmaths not found in PATH (FSL not set up?)" >&2
  exit 1
fi
case "${CHARM_BRAIN_MASK_MODE}" in
  charm|hcp) ;;
  *)
    echo "ERROR: CHARM_BRAIN_MASK_MODE must be charm|hcp (got ${CHARM_BRAIN_MASK_MODE})" >&2
    exit 1
    ;;
esac
case "${CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS}" in
  0|1) ;;
  *)
    echo "ERROR: CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS must be 0|1 (got ${CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS})" >&2
    exit 1
    ;;
esac
case "${CHARM_WRITE_CORTICAL_RIBBON}" in
  0|1) ;;
  *)
    echo "ERROR: CHARM_WRITE_CORTICAL_RIBBON must be 0|1 (got ${CHARM_WRITE_CORTICAL_RIBBON})" >&2
    exit 1
    ;;
esac
case "${CHARM_SKIP_RUN}" in
  0|1) ;;
  *)
    echo "ERROR: CHARM_SKIP_RUN must be 0|1 (got ${CHARM_SKIP_RUN})" >&2
    exit 1
    ;;
esac
if ! [[ "${CHARM_BRAIN_MASK_DILATE_ITERS}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CHARM_BRAIN_MASK_DILATE_ITERS must be an integer >= 0 (got ${CHARM_BRAIN_MASK_DILATE_ITERS})" >&2
  exit 1
fi

# ---- run CHARM ----
echo "Preparing CHARM outputs (logging to ${LOG_FILE})"
cd "${ANAT_DIR}" || exit 1

if [ "${CHARM_SKIP_RUN}" = "1" ]; then
  if [ -d "${ANAT_DIR}/charm" ]; then
    echo "Skipping CHARM run; using existing output folder: ${ANAT_DIR}/charm" >&2
    echo "Skipped CHARM execution on $(date); reused existing anat/charm outputs." > "${LOG_FILE}"
  elif [ -d "${ANAT_DIR}/m2m_${Subject}" ]; then
    echo "Skipping CHARM run; promoting existing folder to anat/charm: ${ANAT_DIR}/m2m_${Subject}" >&2
    rm -rf "${ANAT_DIR}/charm" >/dev/null 2>&1 || true
    mv "${ANAT_DIR}/m2m_${Subject}" "${ANAT_DIR}/charm"
    echo "Skipped CHARM execution on $(date); promoted existing m2m output." > "${LOG_FILE}"
  else
    echo "ERROR: CHARM_SKIP_RUN=1 but neither ${ANAT_DIR}/charm nor ${ANAT_DIR}/m2m_${Subject} exists." >&2
    exit 1
  fi
else
  # remove stale target folder so mv is deterministic
  rm -rf "${ANAT_DIR}/charm" >/dev/null 2>&1 || true

  if [ "${CHARM_REUSE_EXISTING_M2M}" = "1" ] && [ -f "${ANAT_DIR}/m2m_${Subject}/final_tissues.nii.gz" ]; then
    echo "Reusing existing CHARM output folder: ${ANAT_DIR}/m2m_${Subject}" >&2
    echo "Reused existing CHARM output on $(date)" > "${LOG_FILE}"
  else
    # ---- preflight: ensure environment allows local sockets for MPI ----
    if ! python3 - <<'PY'
import socket
try:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    s.listen(1)
    s.close()
except Exception:
    raise SystemExit(1)
PY
    then
      echo "ERROR: environment blocks MPI/socket" >&2
      echo "       CHARM cannot start in this runtime. Re-run in an environment that allows local sockets." >&2
      exit 1
    fi

    # remove any stale output folder before re-running CHARM
    rm -rf "${ANAT_DIR}/m2m_${Subject}" >/dev/null 2>&1 || true

    charm_rc=0
    set +e
    if [ -f "${CHARM_T2}" ]; then
      "${CHARM_BIN}" "${Subject}" "${CHARM_T1}" "${CHARM_T2}" \
        --skipregisterT2 --forceqform --noneck \
        > "${LOG_FILE}" 2>&1
      charm_rc=$?
    else
      "${CHARM_BIN}" "${Subject}" "${CHARM_T1}" \
        --forceqform --noneck \
        > "${LOG_FILE}" 2>&1
      charm_rc=$?
    fi
    set -e

    if [ "${charm_rc}" -ne 0 ] && [ ! -f "${ANAT_DIR}/m2m_${Subject}/final_tissues.nii.gz" ]; then
      echo "ERROR: CHARM exited with code ${charm_rc} and did not produce final_tissues.nii.gz" >&2
      echo "       See ${LOG_FILE} for details." >&2
      exit 1
    fi
    if [ "${charm_rc}" -ne 0 ] && [ -f "${ANAT_DIR}/m2m_${Subject}/final_tissues.nii.gz" ]; then
      echo "WARNING: CHARM exited with code ${charm_rc} but produced outputs; continuing with post-processing." >&2
    fi
  fi

  # Rename m2m_<Subject> --> charm and move log
  if [ -d "${ANAT_DIR}/m2m_${Subject}" ]; then
    mv "${ANAT_DIR}/m2m_${Subject}" "${ANAT_DIR}/charm"
    mv "${LOG_FILE}" "${ANAT_DIR}/charm/charm_log.txt" 2>/dev/null || cp -f "${ANAT_DIR}/charm/charm_log.html" "${ANAT_DIR}/charm/charm_log.txt" 2>/dev/null || true
  else
    echo "ERROR: Expected CHARM output folder not found: ${ANAT_DIR}/m2m_${Subject}" >&2
    echo "       See ${LOG_FILE} for details (if it exists)." >&2
    exit 1
  fi
fi

# Verify required CHARM outputs after either run or reuse/skip path.
if [ ! -f "${CHARM_TISSUE}" ]; then
  echo "ERROR: Missing CHARM tissue file: ${CHARM_TISSUE}" >&2
  exit 1
fi
if [ ! -f "${CHARM_LABELING}" ]; then
  echo "ERROR: Missing CHARM labeling file: ${CHARM_LABELING}" >&2
  exit 1
fi

# ---- aparc+aseg cleanup ----
echo "Cleaning aparc+aseg using CHARM tissues, then CHARM labeling ROI retention"

if [ ! -f "${ASEG_FILE}" ]; then
  echo "ERROR: Missing aseg file: ${ASEG_FILE}" >&2
  exit 1
fi

cd "${T1W_DIR}" || exit 1

# Preserve original aparc+aseg (only create once)
if [ ! -f "aparc+aseg_original.nii.gz" ]; then
  cp "aparc+aseg.nii.gz" "aparc+aseg_original.nii.gz"
fi
if [ ! -f "T1w_acpc_brain_mask.precharm.nii.gz" ] && [ -f "T1w_acpc_brain_mask.nii.gz" ]; then
  cp -f "T1w_acpc_brain_mask.nii.gz" "T1w_acpc_brain_mask.precharm.nii.gz"
fi
if [ ! -f "T1w_acpc_brain.precharm.nii.gz" ] && [ -f "T1w_acpc_brain.nii.gz" ]; then
  cp -f "T1w_acpc_brain.nii.gz" "T1w_acpc_brain.precharm.nii.gz"
fi
if [ ! -f "T1w_acpc_dc_restore_brain.precharm.nii.gz" ] && [ -f "T1w_acpc_dc_restore_brain.nii.gz" ]; then
  cp -f "T1w_acpc_dc_restore_brain.nii.gz" "T1w_acpc_dc_restore_brain.precharm.nii.gz"
fi
if [ ! -f "T2w_acpc_dc_restore_brain.precharm.nii.gz" ] && [ -f "T2w_acpc_dc_restore_brain.nii.gz" ]; then
  cp -f "T2w_acpc_dc_restore_brain.nii.gz" "T2w_acpc_dc_restore_brain.precharm.nii.gz"
fi

orig="aparc+aseg_original.nii.gz"
lab="${CHARM_LABELING}"
tiss="${CHARM_TISSUE}"
target="aparc+aseg.nii.gz"

# Step 1: apply legacy exclusion mask from CHARM final_tissues (before reinsertion).
fslmaths "${tiss}" -thr 1 -uthr 1 -bin tmp_tiss_1.nii.gz
fslmaths "${tiss}" -thr 3 -uthr 3 -bin tmp_tiss_3.nii.gz
fslmaths "${tiss}" -thr 5 -uthr 5 -bin tmp_tiss_5.nii.gz
fslmaths "${tiss}" -thr 7 -uthr 7 -bin tmp_tiss_7.nii.gz
fslmaths "${tiss}" -thr 9 -uthr 9 -bin tmp_tiss_9.nii.gz
fslmaths tmp_tiss_1.nii.gz -add tmp_tiss_3.nii.gz -add tmp_tiss_5.nii.gz -add tmp_tiss_7.nii.gz -add tmp_tiss_9.nii.gz -bin tmp_tiss_cleanup_mask.nii.gz
fslmaths tmp_tiss_cleanup_mask.nii.gz -binv tmp_tiss_keep_mask.nii.gz
fslmaths "${target}" -mas tmp_tiss_keep_mask.nii.gz "${target}"

# Step 2: targeted reinsertion/edits using CHARM labeling.
# pallidum from agreement between original aparc+aseg and CHARM labeling
fslmaths "${orig}" -thr 13 -uthr 13 -bin tmp_orig_13.nii.gz
fslmaths "${lab}"  -thr 13 -uthr 13 -bin tmp_lab_13.nii.gz
fslmaths tmp_orig_13.nii.gz -mul tmp_lab_13.nii.gz -mul 13 tmp_add_13.nii.gz

fslmaths "${orig}" -thr 52 -uthr 52 -bin tmp_orig_52.nii.gz
fslmaths "${lab}"  -thr 52 -uthr 52 -bin tmp_lab_52.nii.gz
fslmaths tmp_orig_52.nii.gz -mul tmp_lab_52.nii.gz -mul 52 tmp_add_52.nii.gz

fslmaths "${target}" -thr 13 -uthr 13 -bin tmp_is_13.nii.gz
fslmaths "${target}" -thr 52 -uthr 52 -bin tmp_is_52.nii.gz
fslmaths tmp_is_13.nii.gz -add tmp_is_52.nii.gz -bin tmp_is_pallidum.nii.gz
fslmaths tmp_is_pallidum.nii.gz -binv tmp_not_pallidum.nii.gz
fslmaths "${target}" -mul tmp_not_pallidum.nii.gz tmp_work.nii.gz
fslmaths tmp_work.nii.gz -add tmp_add_13.nii.gz -add tmp_add_52.nii.gz tmp_work.nii.gz

# brainstem + diencephalon from CHARM labeling
fslmaths "${lab}" -thr 16 -uthr 16 -bin -mul 16 tmp_add_16.nii.gz
fslmaths "${lab}" -thr 28 -uthr 28 -bin -mul 28 tmp_add_28.nii.gz
fslmaths "${lab}" -thr 60 -uthr 60 -bin -mul 60 tmp_add_60.nii.gz

fslmaths tmp_work.nii.gz -thr 16 -uthr 16 -bin tmp_is_16.nii.gz
fslmaths tmp_work.nii.gz -thr 28 -uthr 28 -bin tmp_is_28.nii.gz
fslmaths tmp_work.nii.gz -thr 60 -uthr 60 -bin tmp_is_60.nii.gz
fslmaths tmp_is_16.nii.gz -add tmp_is_28.nii.gz -add tmp_is_60.nii.gz -bin tmp_old_bs_dien.nii.gz
fslmaths tmp_old_bs_dien.nii.gz -binv tmp_not_old_bs_dien.nii.gz
fslmaths tmp_work.nii.gz -mul tmp_not_old_bs_dien.nii.gz tmp_work.nii.gz
fslmaths tmp_work.nii.gz -add tmp_add_16.nii.gz -add tmp_add_28.nii.gz -add tmp_add_60.nii.gz tmp_work.nii.gz

# remove WM coded as 2
fslmaths tmp_work.nii.gz -thr 2 -uthr 2 -bin tmp_is_2.nii.gz
fslmaths tmp_is_2.nii.gz -binv tmp_not_2.nii.gz
fslmaths tmp_work.nii.gz -mul tmp_not_2.nii.gz "${target}"

# Delete temporary/intermediate files
rm -f tmp_orig_13.nii.gz tmp_lab_13.nii.gz tmp_add_13.nii.gz \
  tmp_orig_52.nii.gz tmp_lab_52.nii.gz tmp_add_52.nii.gz \
  tmp_is_13.nii.gz tmp_is_52.nii.gz tmp_is_pallidum.nii.gz tmp_not_pallidum.nii.gz \
  tmp_work.nii.gz tmp_add_16.nii.gz tmp_add_28.nii.gz tmp_add_60.nii.gz \
  tmp_is_16.nii.gz tmp_is_28.nii.gz tmp_is_60.nii.gz tmp_old_bs_dien.nii.gz \
  tmp_not_old_bs_dien.nii.gz tmp_is_2.nii.gz tmp_not_2.nii.gz \
  tmp_tiss_1.nii.gz tmp_tiss_3.nii.gz tmp_tiss_5.nii.gz tmp_tiss_7.nii.gz tmp_tiss_9.nii.gz \
  tmp_tiss_cleanup_mask.nii.gz tmp_tiss_keep_mask.nii.gz

# ---- create binarized brain mask from CHARM labeling (<100 = brain) ----
echo "Creating T1w brain mask from CHARM labeling"

if [ ! -f "${CHARM_LABELING}" ]; then
  echo "ERROR: Missing CHARM labeling file: ${CHARM_LABELING}" >&2
  exit 1
fi

# Anything < 100 is "brain"
CHARM_MASK_TMP="${T1W_DIR}/T1w_acpc_brain_mask.charmtmp.nii.gz"
fslmaths "${CHARM_LABELING}" -uthr 99 -bin "${CHARM_MASK_TMP}"
for ((i = 0; i < CHARM_BRAIN_MASK_DILATE_ITERS; i++)); do
  fslmaths "${CHARM_MASK_TMP}" -dilM "${CHARM_MASK_TMP}"
done
PRECHARM_MASK="${T1W_DIR}/T1w_acpc_brain_mask.precharm.nii.gz"
if [ "${CHARM_BRAIN_MASK_MODE}" = "hcp" ]; then
  if [ ! -f "${PRECHARM_MASK}" ]; then
    echo "WARNING: HCP mode requested but pre-CHARM mask missing; falling back to CHARM mask: ${PRECHARM_MASK}" >&2
    cp -f "${CHARM_MASK_TMP}" "${BRAIN_MASK_OUT}"
  else
    cp -f "${PRECHARM_MASK}" "${BRAIN_MASK_OUT}"
  fi
else
  cp -f "${CHARM_MASK_TMP}" "${BRAIN_MASK_OUT}"
fi

# Ensure mask is binary
fslmaths "${BRAIN_MASK_OUT}" -bin "${BRAIN_MASK_OUT}"
rm -f "${CHARM_MASK_TMP}"

# ---- apply mask to T1w/T2w images ----
echo "Applying brain mask to T1w/T2w volumes"

if [ ! -f "${T1W_ACPC}" ]; then
  echo "ERROR: Missing ${T1W_ACPC}" >&2
  exit 1
fi

if [ ! -f "${T1W_ACPC_DCRESTORE}" ]; then
  echo "ERROR: Missing ${T1W_ACPC_DCRESTORE}" >&2
  exit 1
fi

if [ ! -f "${T2W_ACPC_DCRESTORE}" ]; then
  echo "ERROR: Missing ${T2W_ACPC_DCRESTORE}" >&2
  exit 1
fi

# Apply
fslmaths "${T1W_ACPC}" -mas "${BRAIN_MASK_OUT}" "${T1W_ACPC_BRAIN}"
fslmaths "${T1W_ACPC_DCRESTORE}" -mas "${BRAIN_MASK_OUT}" "${T1W_ACPC_DCRESTORE_BRAIN}"
fslmaths "${T2W_ACPC_DCRESTORE}" -mas "${BRAIN_MASK_OUT}" "${T2W_ACPC_DCRESTORE_BRAIN}"

# ---- create binarized cortical ribbon mask (values 3 and 42) ----
#echo "Creating CorticalRibbon.nii.gz from ribbon.nii.gz (labels 3 and 42)"

RIBBON_FILE="${T1W_DIR}/ribbon.nii.gz"
CORT_RIBBON_OUT="${T1W_DIR}/CorticalRibbon.nii.gz"

if [ ! -f "${RIBBON_FILE}" ]; then
  echo "ERROR: Missing ribbon file: ${RIBBON_FILE}" >&2
  exit 1
fi

# Make binary mask where ribbon == 3 OR ribbon == 42
fslmaths "${RIBBON_FILE}" -thr 3  -uthr 3  -bin "${T1W_DIR}/_ribbon3_tmp.nii.gz"
fslmaths "${RIBBON_FILE}" -thr 42 -uthr 42 -bin "${T1W_DIR}/_ribbon42_tmp.nii.gz"
fslmaths "${T1W_DIR}/_ribbon3_tmp.nii.gz" -add "${T1W_DIR}/_ribbon42_tmp.nii.gz" -bin "${CORT_RIBBON_OUT}"

rm -f "${T1W_DIR}/_ribbon3_tmp.nii.gz" "${T1W_DIR}/_ribbon42_tmp.nii.gz"

# ---- optionally remove HC/Amy + CSF(520) from CorticalRibbon using CHARM labeling ----
if [ "${CHARM_CORTICAL_RIBBON_EXCLUDE_LABELS}" = "1" ]; then
echo "Removing HC/Amy (17,18,53,54) + CSF (520) from CorticalRibbon.nii.gz using CHARM labeling"

CORT_RIBBON_OUT="${T1W_DIR}/CorticalRibbon.nii.gz"
CHARM_LABELING="${ANAT_DIR}/charm/segmentation/labeling.nii.gz"

if [ ! -f "${CORT_RIBBON_OUT}" ]; then
  echo "ERROR: Missing cortical ribbon mask: ${CORT_RIBBON_OUT}" >&2
  exit 1
fi

if [ ! -f "${CHARM_LABELING}" ]; then
  echo "ERROR: Missing CHARM labeling file: ${CHARM_LABELING}" >&2
  exit 1
fi

# Build exclusion mask from CHARM labeling values
# 17/53 = L/R hippocampus, 18/54 = L/R amygdala, 520 = CSF
fslmaths "${CHARM_LABELING}" -thr 17  -uthr 17  -bin "${T1W_DIR}/_hcL_tmp.nii.gz"
fslmaths "${CHARM_LABELING}" -thr 53  -uthr 53  -bin "${T1W_DIR}/_hcR_tmp.nii.gz"
fslmaths "${CHARM_LABELING}" -thr 18  -uthr 18  -bin "${T1W_DIR}/_amyL_tmp.nii.gz"
fslmaths "${CHARM_LABELING}" -thr 54  -uthr 54  -bin "${T1W_DIR}/_amyR_tmp.nii.gz"
fslmaths "${CHARM_LABELING}" -thr 520 -uthr 520 -bin "${T1W_DIR}/_csf520_tmp.nii.gz"

fslmaths "${T1W_DIR}/_hcL_tmp.nii.gz" \
  -add "${T1W_DIR}/_hcR_tmp.nii.gz" \
  -add "${T1W_DIR}/_amyL_tmp.nii.gz" \
  -add "${T1W_DIR}/_amyR_tmp.nii.gz" \
  -add "${T1W_DIR}/_csf520_tmp.nii.gz" \
  -bin "${T1W_DIR}/_exclude_tmp.nii.gz"

# Keep CorticalRibbon voxels that are NOT in the exclusion mask
fslmaths "${T1W_DIR}/_exclude_tmp.nii.gz" -binv "${T1W_DIR}/_keep_tmp.nii.gz"
fslmaths "${CORT_RIBBON_OUT}" -mas "${T1W_DIR}/_keep_tmp.nii.gz" "${CORT_RIBBON_OUT}"

rm -f "${T1W_DIR}/_hcL_tmp.nii.gz" \
      "${T1W_DIR}/_hcR_tmp.nii.gz" \
      "${T1W_DIR}/_amyL_tmp.nii.gz" \
      "${T1W_DIR}/_amyR_tmp.nii.gz" \
      "${T1W_DIR}/_csf520_tmp.nii.gz" \
      "${T1W_DIR}/_exclude_tmp.nii.gz" \
      "${T1W_DIR}/_keep_tmp.nii.gz"
fi

if [ "${CHARM_WRITE_CORTICAL_RIBBON}" = "0" ]; then
  rm -f "${CORT_RIBBON_OUT}"
fi

echo "Done."

#!/usr/bin/env bash
# Config for raw DICOM -> RevisedMe-fMRIPipeline intake import.
#
# Usage:
#   bash bin/mefmri_import_raw.sh <RawDicomDir> <SubjectDir> [ConfigFile] [--session N] [--dry-run]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

IMPORT_PROTOCOL_NAME="cmrr_mb4_te4_rest_default"
IMPORT_DCM2NIIX_BIN="${IMPORT_DCM2NIIX_BIN:-dcm2niix}"

FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"

# Per-import expectations.
IMPORT_EXPECT_REST_RUNS_PER_SESSION=2
IMPORT_EXPECT_ECHOES_PER_RUN=4
IMPORT_EXPECT_SBREF_PER_RUN=1
IMPORT_EXPECT_FMAP_AP_PER_SESSION=2
IMPORT_EXPECT_FMAP_PA_PER_SESSION=2
IMPORT_EXPECT_T1W_MAX_PER_IMPORT=1
IMPORT_EXPECT_T2W_MAX_PER_IMPORT=1

# Session-1 anatomical policy.
IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING=1
IMPORT_T2W_OPTIONAL=1

# Volume-count guardrails.
IMPORT_EXPECT_REST_VOLUMES=650
IMPORT_EXPECT_SBREF_VOLUMES=1
IMPORT_EXPECT_FMAP_VOLUMES=1
IMPORT_EXPECT_ANAT_VOLUMES=1

# File-size guardrails (bytes). These should catch obviously broken outputs.
IMPORT_MIN_BYTES_REST=200000000
IMPORT_MIN_BYTES_SBREF=200000
IMPORT_MIN_BYTES_FMAP=200000
IMPORT_MIN_BYTES_T1W=10000000
IMPORT_MIN_BYTES_T2W=10000000

# Classification rules.
IMPORT_T1W_REGEX='SAG 3D MPRAGE 4 ECHO RMS$'
IMPORT_T2W_REGEX='SAG T2w_0\.8mm$'
IMPORT_REST_REGEX='CMRR-MB4-TE4-2\.5-REST[0-9]+$'
IMPORT_SBREF_REGEX='CMRR-MB4-TE4-2\.5-REST[0-9]+_SBRef$'
IMPORT_FMAP_AP_REGEX='SpinEchoFieldMap_AP_Run[0-9]+$'
IMPORT_FMAP_PA_REGEX='SpinEchoFieldMap_PA_Run[0-9]+$'

# Ignore these series descriptions during classification.
IMPORT_IGNORE_REGEXES=(
  '^AAHead_Scout'
  '_Pha$'
  '_MPR_cor$'
  '_MPR_tra$'
  '^SAG 3D MPRAGE 4 ECHO$'
)

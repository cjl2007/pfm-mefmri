#!/usr/bin/env bash
# Config for raw DICOM -> MultiEchofMRI pipeline intake import.
#
# This is a study-facing template. Update the protocol name, expected counts,
# and regex rules to match your scanner export naming.

CONFIG_DEFAULT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MEDIR="${PFM_MEFMRI_ROOT:-${MEDIR:-$CONFIG_DEFAULT_ROOT}}"

# Copy outside the repository, fill in your acquisition rules, then set to 1.
IMPORT_CONFIG_READY=0
IMPORT_PROTOCOL_NAME=""
IMPORT_DCM2NIIX_BIN="${IMPORT_DCM2NIIX_BIN:-dcm2niix}"

FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"

# Per-import expectations.
IMPORT_EXPECT_TASK_RUNS_PER_SESSION=0
IMPORT_EXPECT_ECHOES_PER_RUN=0
IMPORT_EXPECT_SBREF_PER_RUN=0
IMPORT_EXPECT_FMAP_AP_PER_SESSION=0
IMPORT_EXPECT_FMAP_PA_PER_SESSION=0
IMPORT_EXPECT_T1W_MAX_PER_IMPORT=0
IMPORT_EXPECT_T2W_MAX_PER_IMPORT=0

# Session-1 anatomical policy.
IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING=1
IMPORT_T2W_OPTIONAL=1

# Volume-count guardrails.
IMPORT_EXPECT_TASK_VOLUMES=0
IMPORT_EXPECT_SBREF_VOLUMES=1
IMPORT_EXPECT_FMAP_VOLUMES=1
IMPORT_EXPECT_ANAT_VOLUMES=1

# File-size guardrails (bytes). Set these for your sequence; zero disables each minimum-size check.
IMPORT_MIN_BYTES_TASK=0
IMPORT_MIN_BYTES_SBREF=0
IMPORT_MIN_BYTES_FMAP=0
IMPORT_MIN_BYTES_T1W=0
IMPORT_MIN_BYTES_T2W=0

# Classification rules.
IMPORT_T1W_REGEX='(?!)'
IMPORT_T2W_REGEX='(?!)'
IMPORT_TASK_REGEX='(?!)'
IMPORT_SBREF_REGEX='(?!)'
IMPORT_FMAP_AP_REGEX='(?!)'
IMPORT_FMAP_PA_REGEX='(?!)'

# Ignore these series descriptions during classification.
IMPORT_IGNORE_REGEXES=()

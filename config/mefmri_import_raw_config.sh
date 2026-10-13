#!/usr/bin/env bash
# Compatibility filename; new studies should copy mefmri_import_config.sh.
# Study-specific raw DICOM intake. Copy and edit; see docs/configuration.md.
# Counts below are exact expectations, not universal acquisition defaults.
# Fill in the blank counts and scanner SeriesDescription regexes, then set 1.
IMPORT_CONFIG_READY=0
IMPORT_PROTOCOL_NAME=""
IMPORT_DCM2NIIX_BIN="${IMPORT_DCM2NIIX_BIN:-dcm2niix}"

# Match these to the processing config.
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"

# Required acquisition-specific counts. 0 means none expected, not unchecked.
IMPORT_EXPECT_TASK_RUNS_PER_SESSION=""
IMPORT_EXPECT_ECHOES_PER_RUN=""
IMPORT_EXPECT_SBREF_PER_RUN=""
IMPORT_EXPECT_FMAP_AP_PER_SESSION=""
IMPORT_EXPECT_FMAP_PA_PER_SESSION=""
IMPORT_EXPECT_T1W_MAX_PER_IMPORT=1
IMPORT_EXPECT_T2W_MAX_PER_IMPORT=1
IMPORT_REQUIRE_T1W_IF_SUBJECT_MISSING=1
IMPORT_T2W_OPTIONAL=1

# Task volumes: 0 disables this check; positive values check the expected length.
IMPORT_EXPECT_TASK_VOLUMES=0
IMPORT_EXPECT_SBREF_VOLUMES=1
IMPORT_EXPECT_FMAP_VOLUMES=1
IMPORT_EXPECT_ANAT_VOLUMES=1
IMPORT_ECHO_DIM4_POLICY="abort"  # abort|truncate_to_min; retain abort unless intentional

# Optional file-size minimums in bytes; 0 disables the size check.
IMPORT_MIN_BYTES_TASK=0
IMPORT_MIN_BYTES_SBREF=0
IMPORT_MIN_BYTES_FMAP=0
IMPORT_MIN_BYTES_T1W=0
IMPORT_MIN_BYTES_T2W=0

# Python regexes against scanner SeriesDescription. (?!) matches nothing.
# Replace for acquired types; retain (?!) for absent types. Never use an empty
# regex: it matches every series. Example only: IMPORT_TASK_REGEX='^MY_REST$'
IMPORT_T1W_REGEX='(?!)'
IMPORT_T2W_REGEX='(?!)'
IMPORT_TASK_REGEX='(?!)'
IMPORT_SBREF_REGEX='(?!)'
IMPORT_FMAP_AP_REGEX='(?!)'
IMPORT_FMAP_PA_REGEX='(?!)'
IMPORT_IGNORE_REGEXES=()  # Add only series you intentionally exclude.

if [[ "$IMPORT_CONFIG_READY" != 1 ]]; then
  echo 'Edit the raw-import config for your acquisition, then set IMPORT_CONFIG_READY=1. See docs/configuration.md.' >&2
  return 2
fi
for import_count_key in IMPORT_EXPECT_TASK_RUNS_PER_SESSION IMPORT_EXPECT_ECHOES_PER_RUN \
  IMPORT_EXPECT_SBREF_PER_RUN IMPORT_EXPECT_FMAP_AP_PER_SESSION IMPORT_EXPECT_FMAP_PA_PER_SESSION; do
  if [[ ! "${!import_count_key}" =~ ^[0-9]+$ ]]; then
    echo "Set $import_count_key to the expected acquisition count (0 means none)." >&2
    return 2
  fi
done
unset import_count_key

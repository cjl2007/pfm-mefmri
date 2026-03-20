#!/usr/bin/env bash
# Container-oriented preset.
# Usage:
#   bash bin/mefmri_pipeline_container.sh <SubjectDir> config/mefmri_wrapper_config_container.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mefmri_wrapper_config.sh"

# In container images we install tedana directly into PATH.
TEDANA_ACTIVATE_MODE="direct"
TEDANA_ENV=""

# Typical in-container location when SimNIBS is installed.
# Leave empty to pass through default CHARM behavior from module arguments.
CHARM_BIN="${CHARM_BIN:-/opt/simnibs/bin/charm}"

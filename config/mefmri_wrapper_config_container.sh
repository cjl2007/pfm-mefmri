#!/usr/bin/env bash
# Start from the portable study template, then configure tools on container PATH.
# The container launcher supplies the repository root; use its bundled template.
source "${PFM_MEFMRI_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}/config/mefmri_wrapper_config.sh"
TEDANA_ACTIVATE_MODE="direct"
AROMA_ACTIVATE_MODE="direct"
WARPKIT_ACTIVATE_MODE="direct"
# Install optional branch dependencies into the image before enabling them.

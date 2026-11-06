# TANS-only CONFIGURATION FRAGMENT, not a complete pfm-mefmri config or launcher.
# Copy into an existing study config, or source a study-owned copy at its end.
# Edit all installation/study paths. See docs/tans-integration.md here and the
# full configuration reference in the separate pfm-tans repository.

TANS_ENABLE=1
TANS_CODE_ROOT="/path/to/pfm-tans/pipeline"
TANS_CONFIG_FILE="/path/to/MyStudy/config/tans.yaml"
TANS_PYTHON="/path/to/SimNIBS-4.5/simnibs_env/bin/python"
TANS_OUTPUT_DIRNAME="tans"  # One folder name under SubjectDir; no slash.
TANS_PFM_DIR=""            # Default: <SubjectDir>/pfm/RidgePFM; legacy locations remain readable.
TANS_WORKERS=6

TANS_TARGET_NETWORK="" # Set the network label for your project.
TANS_TARGET_SEARCH_SPACE="glasser_predefined_plus_subject_specific_mfg"
TANS_COARSE_GRID_DENSITY=0.85
TANS_INCLUDE_SM_HAND_CALIBRATION=1
TANS_ROBUSTNESS_ENABLE=1    # 0 selects nominal fine winner; robustness not evaluated.

# Example hardware: verify against your physical setup, including sibling STL.
TANS_SM_HAND_COIL_MODEL="/path/to/SimNIBS-4.5/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_C-B60.ccd"
TANS_TREATMENT_COIL_MODEL="/path/to/SimNIBS-4.5/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B65.ccd"
TANS_CALIBRATION_REGISTRY=""  # Default: bundled registry.
TANS_SM_HAND_DIDT_CALIBRATION="magventure_x100_c-b60_biphasic_normal"
TANS_TREATMENT_DIDT_CALIBRATION="magventure_x100_cool-b65_biphasic_normal"
TANS_ALLOW_SIMILAR_COIL_CALIBRATION=0
TANS_CALIBRATION_SUBSTITUTION_JUSTIFICATION=""
TANS_MAXIMUM_TREATMENT_RMT_RATIO=1.20  # Explicitly empty disables the upward bound.

# Supported magnitude-equivalent coils only; review before changing hardware.
TANS_HANDLE_ORIENTATION_PREFERENCE="posterior_inferior"
TANS_COIL_180_DEGREE_MAGNITUDE_EQUIVALENT=1

TANS_START_FROM="headmodel"
TANS_STOP_AFTER="reports"
TANS_RESUME=1

# Set these in the COMPLETE pfm-mefmri config (not a TANS YAML):
# STOP_AFTER_MODULE="tans"
# For existing, complete anatomy/PFM only: START_FROM_MODULE="tans"
# For the full upstream chain: retain the appropriate upstream start stage.

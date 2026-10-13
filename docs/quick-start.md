# Quick Start

This page gives the shortest path from a prepared subject directory to pipeline
outputs. See [Pipeline Overview](pipeline-overview.md) for the full module map.

## 1. Install Dependencies

Install the tools listed in `SOFTWARE_DEPENDENCIES.txt`, then make sure the
configured environments expose the expected commands. Common requirements are:

- HCP Pipelines and Connectome Workbench (`wb_command`).
- FSL commands such as `fslmaths`, `flirt`, `applywarp`, and `topup`.
- FreeSurfer with a readable `license.txt`.
- An MSM build if using `HCP_REGNAME="MSMSulc"`.
- SimNIBS CHARM if using CHARM-derived anatomical masks.
- tedana/ME-ICA, ICA-AROMA, MEDIC/warpkit, and `pfm-nsi` when those branches
  are enabled.

## 2. Configure The Pipeline

Follow [configuration setup](configuration.md) to copy the processing and site
environment templates into your study folder. Set your software paths and pass
that study config explicitly when launching; do not run the config as a script.
The [dependency guide](dependencies.md) separates core and optional tools.

Review the path and environment settings near the top first:

```bash
EnvironmentScript=/path/to/MyStudy/config/environment.sh
CHARM_BIN=/path/to/charm
FS_LICENSE=/path/to/license.txt
```

Then set the main routing options:

```bash
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"
PROCESSING_MODE="auto"
DISTORTION_CORRECTION_MODE="medic"
MULTI_ECHO_DENOISE_METHOD="meica"
SINGLE_ECHO_DENOISE_METHOD="acompcor"
PFM_ENABLE=1  # optional; leave 0 for preprocessing/QC/NSI only
PFM_STRATEGY="ridgepfm"
```

For a fresh run, use:

```bash
START_FROM_MODULE="validate"
STOP_AFTER_MODULE=""
```

For debugging or reruns, `START_FROM_MODULE` and `STOP_AFTER_MODULE` accept:

```text
validate | anat_hcp | anat_charm | fieldmaps | coreg | headmotion |
denoise | mgtr | vol2surf | concat | nsi | pfm | pfm_update | tans
```

When a study requests MEDIC but a run uses the automatic TOPUP fallback,
restarting at `coreg` or `headmotion` resolves the mode from existing outputs.
The runner retains MEDIC when every selected run has its expected nonempty
displacement map. If none has a MEDIC map, it can reuse the saved TOPUP outputs:
the measured-fieldmap marker, acquisition parameters, and all three average
fieldmap images. This does not require repeating raw-phase preflight or TOPUP
estimation. Partial MEDIC outputs or incomplete TOPUP outputs stop the restart
with an actionable error; direct-B0 and no-fieldmap placeholders are not treated
as TOPUP. Run metadata records both the requested and resolved modes.

To change distortion-correction methods, restart at `fieldmaps`. After changing
coregistration, rerun headmotion and downstream analyses so they use the new
transforms. Stop any earlier process using the same subject outputs before
starting a replacement run.

## 3. Prepare Input Layout

The default rest layout is:

```text
<SubjectDir>/
  anat/
    unprocessed/
      T1w/
      T2w/
  func/
    unprocessed/
      rest/
        field_maps/
          AP_S1_R1.nii.gz
          AP_S1_R1.json
          PA_S1_R1.nii.gz
          PA_S1_R1.json
        session_1/
          run_1/
            Rest_S1_R1_E1.nii.gz
            Rest_S1_R1_E1.json
            Rest_S1_R1_E2.nii.gz
            Rest_S1_R1_E2.json
            Rest_S1_R1_E3.nii.gz
            Rest_S1_R1_E3.json
            SBref_S1_R1_E1.nii.gz
            SBref_S1_R1_E1.json
```

Task-local field maps are read from:

```text
func/unprocessed/<FUNC_DIRNAME>/field_maps/
```

Processed field maps are written to:

```text
func/<FUNC_DIRNAME>/field_maps/
```

This prevents one task or acquisition from silently reusing another task's field
maps.

## 4. Run One Subject

```bash
bash bin/mefmri_pipeline.sh /path/to/MyStudy/SUB001 /path/to/MyStudy/config/mefmri.sh
```

Start with `STOP_AFTER_MODULE="validate"` and review validation output before
setting it to `""` for a complete run. The subject directory is the subject folder
itself, not the study root. Omitting the config argument uses the repository template.

## 5. Run All Task Folders

To process every task folder under `func/unprocessed/`:

```bash
FUNC_DIRNAME="All"
```

Anatomy runs once, then each discovered task folder gets its own functional
field maps, transforms, QC, concatenation, NSI, and PFM outputs.

## 6. First Outputs To Review

After a full run, start with:

```text
func/<FUNC_DIRNAME>/QC/Validation/
func/<FUNC_DIRNAME>/QC/Coreg/functional_coregistration_qc.html
func/<FUNC_DIRNAME>/QC/Fieldmaps/distortion_correction_qc.html
func/<FUNC_DIRNAME>/QC/MotionGrayPlots/motion_grayplots_qc.html
func/<FUNC_DIRNAME>/QC/NSI/
func/<FUNC_DIRNAME>/PFM/
```

Exact availability depends on enabled modules and the selected denoising/PFM
strategy.

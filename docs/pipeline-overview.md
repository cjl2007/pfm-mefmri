# Pipeline Overview

PFM-ME fMRI is organized as a staged shell workflow with Python/MATLAB helper
programs for denoising, QC, CIFTI manipulation, NSI, and PFM. The main entry
point is:

```bash
bash bin/mefmri_pipeline.sh <SubjectDir> [ConfigFile]
```

The default config is:

```bash
config/mefmri_wrapper_config.sh
```

## Stage Order

The pipeline stage order is:

```text
validate
anat_hcp
anat_charm
fieldmaps
coreg
headmotion
denoise
mgtr
vol2surf
concat
nsi
pfm
pfm_update
tans
```

Use `START_FROM_MODULE` and `STOP_AFTER_MODULE` to resume or inspect a subset.
The legacy `meica` tag is accepted as an alias for `denoise`.

When `CRAWLING_SEED_FC_ENABLE=1`, the CrawlingSeedFC movie runs as an optional
concat-stage QC sub-step after concat and before NSI.

## Anatomy

The anatomical branch produces HCP-style anatomy and surface outputs used by
functional registration, volume-to-surface mapping, CIFTI generation, PFM, and
Workbench review.

Relevant controls:

```bash
HCP_REGNAME="MSMSulc"       # MSMSulc|FS
CHARM_BRAIN_MASK_MODE="charm"
CHARM_WRITE_CORTICAL_RIBBON=1
```

The CHARM branch can generate a CHARM-derived brain mask and a cortical ribbon
mask. The cortical ribbon is used by MGTR and volume-to-surface mapping when
enabled.

## Functional Routing

Functional naming is controlled by:

```bash
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"
FUNC_XFMS_DIRNAME=""
```

`FUNC_DIRNAME` selects the task folder under `func/unprocessed/`. Set it to
`All` to process all discovered task folders.

`FUNC_XFMS_DIRNAME` normally follows `FUNC_DIRNAME`. Set it only when a task
needs to share an existing transform namespace.

## Importing Inputs

Use `bin/mefmri_import.sh` as the common front door for raw scanner exports and
BIDS datasets. It dispatches to the raw importer when DICOM-like files are
present and to the BIDS importer when BIDS-style `sub-*` NIfTIs are present:

```bash
bash bin/mefmri_import.sh /path/to/raw_dicom_export /path/to/study/ME001 config/mefmri_import_raw_config.sh --session 1 --nordic
bash bin/mefmri_import.sh /path/to/bids 06 /path/to/study/ME06 --task rest --mode symlink --overwrite
```

Use `--input-type raw` or `--input-type bids` when auto-detection is ambiguous.
Raw imports can run/stage NORDIC with `--nordic`; BIDS imports preserve NORDIC
metadata already present in sidecars and write run-level `NORDIC_DENOISING.txt`
markers.

## Distortion Correction

Choose the branch with:

```bash
DISTORTION_CORRECTION_MODE="topup"     # topup|direct_b0|medic|none
```

Supported modes:

- `topup`: paired AP/PA spin-echo EPI field maps.
- `direct_b0`: gradient-echo/direct B0 field maps, including legacy and
  BIDS-style phase/magnitude patterns.
- `medic`: MEDIC/warpkit correction from per-echo phase companions. If MEDIC
  inputs are incomplete and AP/PA field maps are available, the runner can fall
  back to TOPUP.
- `none`: writes zero-unwarp placeholders and skips susceptibility correction.
  Use explicitly only when needed because registration quality can suffer.

## Denoising Branches

The pipeline can process multi-echo and single-echo inputs.

Main controls:

```bash
PROCESSING_MODE="auto"              # auto|multi_echo|single_echo
MULTI_ECHO_DENOISE_METHOD="meica"   # meica|acompcor|aroma
SINGLE_ECHO_DENOISE_METHOD="acompcor"
SINGLE_ECHO_ECHO_INDEX=1
```

Multi-echo options:

- `meica`: tedana/ME-ICA with optional NSI/spatial reclassification.
- `acompcor`: optimally combined multi-echo data followed by aCompCor.
- `aroma`: optimally combined multi-echo data followed by ICA-AROMA.

Single-echo options:

- `acompcor`: E1 or configured source echo followed by aCompCor.
- `aroma`: E1 or configured source echo followed by ICA-AROMA.

## MGTR, Vol2Surf, And CIFTI

MGTR can regress the cortical ribbon mean gray-timecourse after denoising:

```bash
MGTR_ENABLE="auto"      # auto|0|1
```

In `auto` mode, MGTR is skipped after aCompCor because the aCompCor nuisance
model already includes the cortical-ribbon mean signal and derivative.

Volume-to-surface mapping can use the CHARM/HCP cortical ribbon:

```bash
VOL2SURF_USE_CORTICAL_RIBBON_MASK=1
VOL2SURF_USE_GOOD_VOXELS_MASK=0
```

Concatenation combines run-level CIFTIs and writes censored subject-level
dtseries files:

```bash
CONCAT_ENABLE=1
CONCAT_CENSOR_BY_FD=1
CONCAT_FD_THRESHOLD=0.3
```

## Downstream Modules

Downstream analysis modules are independent switches:

```bash
CRAWLING_SEED_FC_ENABLE=1
NSI_ENABLE=1
PFM_ENABLE=1
PFM_STRATEGY="ridgepfm"      # ridgepfm|infomap|all
```

The `pfm_update` stage applies manual Infomap edits without rerunning heavy
community detection.

## Optional TANS Stage

With `TANS_ENABLE=1`, the wrapper calls the Python workflow in `TANS_CODE_ROOT`
after PFM. TANS builds its separate head model from native `anat/T1w/T1w.nii.gz`
and uses `anat/T2w/T2w.nii.gz` when available; without T2w it falls back to a
T1w-only CHARM model. It does not reuse the ACPC `anat/charm` head model.

```bash
TANS_CODE_ROOT="/opt/pfm-tans/pipeline"
TANS_CONFIG_FILE="/data/MyStudy/config/tans.yaml"
TANS_OUTPUT_DIRNAME="tans"
```

The output is always `<SubjectDir>/<TANS_OUTPUT_DIRNAME>`, with an unset or empty
folder name defaulting to `tans`. For example, `TANS_OUTPUT_DIRNAME="tans_review"`
selects `<SubjectDir>/tans_review` for each subject without a subject-specific
config path. Specify one folder name, not an absolute path or nested path.

The old `TANS_OUTPUT_DIR` setting is rejected when nonempty: replace it with
`TANS_OUTPUT_DIRNAME` in existing configs. Changing this setting does not move
or rename existing outputs.

For a subject with completed PFM, set `START_FROM_MODULE="tans"` and
`STOP_AFTER_MODULE="tans"`. `TANS_START_FROM="headmodel"` controls where the
separate TANS workflow begins, not where the ME-fMRI wrapper resumes.

Use the [complete TANS integration guide](tans-integration.md) and portable
shell/YAML examples. Select maps, cortical search region, physical coils, and
calibration identities explicitly for the study.

## Provenance

When enabled:

```bash
RUN_CONFIG_SNAPSHOT=1
```

the pipeline writes a run metadata snapshot under:

```text
func/<FUNC_DIRNAME>/QC/RunMetadata/pipeline_run_<timestamp>.txt
```

This captures effective settings without renaming the canonical outputs.

# Multi-Echo fMRI Pipeline

This repository contains a modular, config-driven fMRI preprocessing and
precision functional mapping pipeline. It is designed for dense longitudinal
sampling of individual participants, where many runs and sessions are combined
to estimate stable subject-specific functional network topographies.

The workflow supports multi-echo and single-echo acquisitions, task-local
distortion correction, denoising, CIFTI generation, concatenation, and precision
functional mapping. Outputs are organized so that an individual participant can
be inspected across sessions, runs, QC modules, and network-mapping products
without losing run-level provenance.

## Release Setup (Required)

Before first run on a new machine:

1. Install the [core and selected optional dependencies](dependencies.md).
2. Follow [configuration setup](configuration.md) to copy and edit your study's
   processing, environment, and (for raw DICOMs) import configurations.
3. If CHARM is not on `PATH`, set `CHARM_BIN` to your local CHARM binary path.
4. Set `FS_LICENSE` or `FS_LICENSE_FILE` to a valid FreeSurfer `license.txt`.

## Config layout

The main config is organized in this order:

1. Core paths and environment settings
2. Common user-facing pipeline and QC options
3. Module-specific settings
4. Advanced and rarely changed settings

This keeps frequent settings near the top and expert knobs near the bottom.

## Running

```bash
bash bin/mefmri_pipeline.sh <SubjectDir> [ConfigFile]
```

Example:

```bash
bash bin/mefmri_pipeline.sh /path/to/MyStudy/SUB001 /path/to/MyStudy/config/mefmri.sh
```

## Importing Inputs

Use the unified importer for either raw scanner exports/DICOMs or BIDS NIfTI
datasets:

```bash
bash bin/mefmri_import.sh <InputDir> <...args...>
```

The importer auto-detects raw DICOM inputs from DICOM-like files (`*.dcm` or
`DICOMDIR`) and BIDS inputs from `sub-*` NIfTIs. You can override detection with
`--input-type raw` or `--input-type bids`.

Examples:

```bash
# Raw scanner export: first copy/edit config/mefmri_import_config.sh as import.sh.
bash bin/mefmri_import.sh /path/to/raw_dicom_export /path/to/study/ME001 /path/to/MyStudy/config/import.sh --session 1

# BIDS root.
bash bin/mefmri_import.sh /path/to/bids 06 /path/to/study/ME06 --task rest --mode symlink --overwrite
```

The backend-specific scripts remain available for advanced/debug use:
`bin/mefmri_import_raw.sh` and `bin/mefmri_import_bids.sh`.

## Expected Input Layout

The default task is resting state:

```bash
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"
```

This is only a convention. Any task folder under `func/unprocessed/` can be
processed independently with the same anatomical outputs by setting
`FUNC_DIRNAME` and `FUNC_FILE_PREFIX` in the config. For example, a task can use
`FUNC_DIRNAME="<TaskName>"` and `FUNC_FILE_PREFIX="<TaskPrefix>"`.

Functional inputs and field maps are task-local:

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
            SBref_S1_R1_E1.nii.gz
            SBref_S1_R1_E1.json
      <TaskName>/
        field_maps/
          AP_S1_R1.nii.gz
          PA_S1_R1.nii.gz
        session_1/
          run_1/
            <TaskPrefix>_S1_R1_E1.nii.gz
            <TaskPrefix>_S1_R1_E1.json
```

The fieldmap stage reads from
`func/unprocessed/<FUNC_DIRNAME>/field_maps/` and writes processed field maps to
`func/<FUNC_DIRNAME>/field_maps/`. This prevents field maps acquired for one
task or acquisition from being reused for another task by accident.

Distortion correction is selected with `DISTORTION_CORRECTION_MODE`:

- `topup`: reverse phase-encoding spin-echo EPI field maps, expected as paired
  `AP_S*_R*.nii.gz` and `PA_S*_R*.nii.gz` files with JSON metadata.
- `direct_b0`: gradient-echo/direct B0 field mapping, expected as paired
  phase/phasediff and magnitude files. Legacy names such as
  `FieldMap_S1_R1.nii.gz` plus `Mag_S1_R1.nii.gz` are supported, as are
  BIDS-style `*_phasediff.nii.gz` plus `*_magnitude1.nii.gz` patterns. Set
  `FM_DIRECT_B0_DELTA_TE_MS` for the `fsl_prepare_fieldmap` delta TE.
- `medic`: MEDIC/warpkit-derived correction path. Requires per-echo phase
  companions; if they are absent and AP/PA field maps are present, the runner
  falls back to `topup`.
- `none`: skip susceptibility correction and write zero-unwarp placeholders.
  Use this only explicitly; it can degrade functional-to-anatomical
  registration.

Set `FUNC_DIRNAME="All"` to process every task folder discovered under
`func/unprocessed/`. Anatomy runs once, then each task is processed separately
with its own functional runs, field maps, transforms, QC, concatenation, and PFM
outputs.

For legacy datasets with shared field maps, you can override:

```bash
FM_RAW_DIR_REL="func/unprocessed/field_maps"
FM_OUT_DIR_REL="func/field_maps"
```

## FreeSurfer License

Set either `FS_LICENSE` or `FS_LICENSE_FILE` to a readable `license.txt` path.

- `anat_hcp` now performs a fail-fast license preflight.

## Tool Path Assumptions (FSL/Workbench/MSM)

The modules call commands like `fslmaths` and `wb_command` by name. This is fine as long as environment setup exposes those binaries on `PATH` (typically via `SetUpHCPPipeline.sh`).

For HCP PostFreeSurfer `--regname=MSMSulc` (the pipeline default), you also need an HCP-compatible MSM (Multimodal Surface Matching) build that supports `--regoption`. Configure this via `MSMBINDIR` (directory containing `msm`, or the msm binary itself). If MSM is not available, set `HCP_REGNAME=FS` in your config to skip MSMSulc registration.

## Provenance

Canonical output filenames are preserved. Provenance is written through logs and run metadata snapshots.

When `RUN_CONFIG_SNAPSHOT=1`, each run writes:

- `func/<FUNC_DIRNAME>/QC/RunMetadata/pipeline_run_<timestamp>.txt`

This captures key settings without renaming canonical outputs.

## PFM Strategies

PFM supports two strategies through `PFM_STRATEGY`.

`infomap` is the community-detection workflow previously used by our group,
including Lynch et al. 2024 Nature. The current Python path preserves that basic
approach, then optionally maps subject-specific Infomap communities onto
canonical network IDs with `PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1`. This mapping
is a post-processing step: it does not force Infomap to return a fixed number of
communities, and multiple communities may receive the same canonical network
label. The labeler combines each community's functional similarity to the
network priors with its spatial overlap, writes confidence maps/tables, flags
low-confidence assignments for review, supports manual correction tables, and
builds mode/probability consensus outputs across graph densities.

`ridgepfm` is a newer, faster PFM strategy. It estimates network assignments
with a ridge-regularized fusion of subject functional connectivity evidence and
spatial/network priors. The prior-guided formulation is intended to reduce the
runtime and improve stability on difficult or lower-quality datasets while still
allowing the subject data to drive deviations from the priors.
For inputs without GSR/MGTR-style global signal control, set
`PFM_RF_FC_DEMEAN=1` to demean each FC fingerprint before RidgePFM estimation.

Both strategies can optionally feed the areal parcellation step
(`PFM_AREAL_ENABLE=1`), which sub-parcellates the network-level output into
smaller areal parcels.

RidgePFM also has a disabled-by-default experimental multi-network layer.
Set `PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE=1` with areal parcellation to test
fixed parcels across blocked random temporal thirds, calibrate added-network
held-out gains with a parcelwise candidate-identity null and cortex/subcortex
FDR, and emit both whole-parcel and localized two-/three-network stripe maps.
The stricter max-statistic result is retained as a high-confidence tier. See
`docs/pfm-outputs.md` for interpretation and controls. Use
`PFM_RF_MULTINETWORK_STRIPE_PRESET="strict|balanced|loose|custom"` as the
single user-facing tightness control; `balanced` is the recommended
experimental starting point. Preliminary behavior is encouraging, but the
feature remains opt-in while it is tested across more participants and
acquisition conditions.

PFM writes strategy-specific outputs:

```text
func/<FUNC_DIRNAME>/PFM/
  RidgePFM/
  Infomap/
```

Infomap uses the Python engine in `lib/pfm_infomap.py`. If `infomap` is not on
`PATH`, set `PFM_INFOMAP_BINARY=/full/path/to/infomap`. For bring-up,
`PFM_INFOMAP_DRY_RUN=1` validates argument wiring without running the heavy
community-detection computation.

Manual Infomap edits are applied with the separate `pfm_update` stage. Enable
`PFM_INFOMAP_UPDATE_ENABLE=1` and run from `START_FROM_MODULE=pfm_update`; the
module scans the per-density manual correction tables, applies edits without
rerunning Infomap, and writes adjusted density and consensus labels.

With areal parcellation enabled, the final areal output is named from the
network-label prefix, for example `InfomapNetworkLabels+ArealParcellation.dlabel.nii`
or `RidgePFM_VTX+ArealParcellation.dlabel.nii`.

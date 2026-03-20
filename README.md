# pfm-mefmri

This repository contains a modular, config-driven multi-echo fMRI preprocessing pipeline for subject directories organized in the expected local project layout.

## Release Defaults

The main pipeline config is `config/mefmri_wrapper_config.sh`.

## Setup

Before first run on a new machine:

1. Clone the repository and initialize submodules:

```bash
git clone https://github.com/cjl2007/pfm-mefmri.git
cd pfm-mefmri
git submodule update --init --recursive
```

2. Install the software listed in `SOFTWARE_DEPENDENCIES.txt`.
3. Ensure FSL, FreeSurfer, Connectome Workbench, GNU Parallel, and Python are available on `PATH`.
4. Review `config/mefmri_wrapper_config.sh` and fill in any machine-specific settings.
5. Set `FS_LICENSE` or `FS_LICENSE_FILE` to a readable FreeSurfer `license.txt`.
6. If CHARM is not on `PATH`, set `CHARM_BIN` in the config.

## pfm-nsi

This release includes `pfm-nsi` as a bundled git submodule at `lib/pfm-nsi`.

Default behavior:

- `NSI_USE_EXTERNAL_CLI=1`
- `NSI_EXTERNAL_ROOT="$MEDIR/lib/pfm-nsi"`

This means a standard clone plus `git submodule update --init --recursive` is enough to use the bundled release-pinned `pfm-nsi`.

If you want to use a different local `pfm-nsi` checkout instead of the bundled one, edit `config/mefmri_wrapper_config.sh` and set:

```bash
NSI_EXTERNAL_ROOT="/full/path/to/your/pfm-nsi"
```

So yes, users can override the bundled release version with their own newer local copy whenever they want.

## Expected Subject Layout

Before preprocessing, each subject directory should be organized as follows:

```text
ME001/
  anat/
    unprocessed/
      T1w/
        T1w_1.nii.gz
        T1w_1.json
        T1w_2.nii.gz
        T1w_2.json
        ...
      T2w/
        T2w_1.nii.gz
        T2w_1.json
        T2w_2.nii.gz
        T2w_2.json
        ...
  func/
    unprocessed/
      rest/
        session_1/
          run_1/
            Rest_S1_R1_E1.nii.gz
            Rest_S1_R1_E1.json
            Rest_S1_R1_E2.nii.gz
            Rest_S1_R1_E2.json
            Rest_S1_R1_E3.nii.gz
            Rest_S1_R1_E3.json
            Rest_S1_R1_E4.nii.gz
            Rest_S1_R1_E4.json
          run_2/
            Rest_S1_R2_E1.nii.gz
            Rest_S1_R2_E1.json
            Rest_S1_R2_E2.nii.gz
            Rest_S1_R2_E2.json
            Rest_S1_R2_E3.nii.gz
            Rest_S1_R2_E3.json
            Rest_S1_R2_E4.nii.gz
            Rest_S1_R2_E4.json
      field_maps/
        AP_S1_R1.nii.gz
        AP_S1_R1.json
        AP_S1_R2.nii.gz
        AP_S1_R2.json
        PA_S1_R1.nii.gz
        PA_S1_R1.json
        PA_S1_R2.nii.gz
        PA_S1_R2.json
```
## Notes

At least one T1-weighted image is required: T1w_1.nii.gz with its matching T1w_1.json.

Additional T1w repeats are supported (T1w_2, T1w_3, etc.) and will be averaged before downstream structural processing and FreeSurfer.

T2-weighted images are optional. If provided, each image should also have a matching JSON sidecar.

Additional T2w repeats are supported and will be averaged in the same way.

Functional runs should be organized by session and run, with one NIfTI file and one matching JSON sidecar per echo.

The number of echoes is flexible. The example above shows 4 echoes, but runs may contain 3, 5, or however many echoes are present in the acquisition.

Field maps should be placed in func/unprocessed/field_maps/, with a matching JSON sidecar for each NIfTI file.

AP and PA are shown for illustrative purposes only. Other phase-encoding direction pairs (for example, RL / LR) are also acceptable, provided they are named consistently and specified correctly in the pipeline configuration.

In general, every input NIfTI file is expected to have a corresponding .json sidecar containing the required acquisition metadata.

## Importers

Two helper entrypoints are provided to build the expected subject layout.

### Import Raw Scanner Exports

Use `bin/mefmri_import_raw.sh` when your input is a raw DICOM export folder.

The default raw import template is `config/mefmri_import_raw_config.sh`. Update its protocol name, expected counts, and regex rules to match your study before use.

Example:

```bash
bash bin/mefmri_import_raw.sh \
  /path/to/raw_dicom_export \
  /path/to/study/ME001 \
  config/mefmri_import_raw_config.sh \
  --session 1
```

Dry-run example:

```bash
bash bin/mefmri_import_raw.sh \
  /path/to/raw_dicom_export \
  /path/to/study/ME001 \
  config/mefmri_import_raw_config.sh \
  --session 1 \
  --dry-run
```

### Import From BIDS

Use `bin/mefmri_import_bids.sh` when your source data are already in BIDS.

This importer maps BIDS inputs into the pipeline's expected local layout and can either symlink or copy files.

Example:

```bash
bash bin/mefmri_import_bids.sh \
  /path/to/bids \
  06 \
  /path/to/study/ME06 \
  --task rest \
  --mode symlink \
  --overwrite
```

## Running The Pipeline

Default invocation:

```bash
bash bin/mefmri_pipeline.sh /path/to/study/ME001
```

Explicit config:

```bash
bash bin/mefmri_pipeline.sh \
  /path/to/study/ME001 \
  config/mefmri_wrapper_config.sh
```

The pipeline default config path is already `config/mefmri_wrapper_config.sh`, so the second argument is only needed when testing an alternate config.

## Provenance And QA

Canonical output filenames are preserved. Provenance is recorded through logs and run metadata snapshots.

When `RUN_CONFIG_SNAPSHOT=1`, each run writes:

- `func/qa/RunMetadata/pipeline_run_<timestamp>.txt`

When `NSI_USE_EXTERNAL_CLI=1`, NSI outputs are written under:

- `func/qa/NSI/`

## Author

Chuck Lynch  
cjl2007@med.cornell.edu

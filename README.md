# pfm-mefmri

This repository contains a modular, config-driven multi-echo fMRI preprocessing pipeline for subject directories organized in the expected local project layout.

## Release Defaults

The main pipeline config is `config/mefmri_wrapper_config.sh`.

## Setup

Before first run on a new machine:

1. Clone the repository and initialize submodules:

```bash
git clone <your-repo-url> pfm-mefmri
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

The pipeline expects each subject directory to look like this before preprocessing:

```text
ME001/
  anat/
    unprocessed/
      T1w/
        T1w_1.nii.gz
      T2w/
        T2w_1.nii.gz
  func/
    unprocessed/
      rest/
        session_1/
          run_1/
            Rest_S1_R1_E1.nii.gz
            Rest_S1_R1_E2.nii.gz
            Rest_S1_R1_E3.nii.gz
            Rest_S1_R1_E4.nii.gz
          run_2/
            Rest_S1_R2_E1.nii.gz
            Rest_S1_R2_E2.nii.gz
            Rest_S1_R2_E3.nii.gz
            Rest_S1_R2_E4.nii.gz
      field_maps/
        AP_S1_R1.nii.gz
        AP_S1_R2.nii.gz
        PA_S1_R1.nii.gz
        PA_S1_R2.nii.gz
```

The pipeline writes preprocessing outputs under the same subject directory, including run metadata and logs in `func/qa/`.

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

# Set up a study configuration

The templates contain processing defaults, not a complete acquisition protocol.
Keep study-specific copies outside the repository. All `/path/to/...` values below
are placeholders to replace with your own paths.

## 1. Get the pipeline and install dependencies

```bash
git clone --recurse-submodules https://github.com/cjl2007/pfm-mefmri.git
cd pfm-mefmri
```

For an existing checkout, run `git submodule update --init --recursive`.
This retrieves the pinned pfm-nsi dependency. It does not install the external
neuroimaging tools. A GitHub ZIP alone does not populate the NSI submodule.

Use the [dependency guide](dependencies.md) to install the tools for your chosen
modules. The [software inventory](../SOFTWARE_DEPENDENCIES.txt) records development
versions and command checks. Do not assume a newer major version gives identical
scientific outputs.

## 2. Copy the configuration files

From the repository root:

```bash
mkdir -p /path/to/MyStudy/config
cp config/mefmri_wrapper_config.sh /path/to/MyStudy/config/mefmri.sh
cp config/mefmri_environment.example.sh /path/to/MyStudy/config/environment.sh
# Only for raw DICOM import:
cp config/mefmri_import_config.sh /path/to/MyStudy/config/import.sh
```

| Study file | Purpose |
| --- | --- |
| `environment.sh` | Machine/site software paths and environment setup. |
| `mefmri.sh` | Processing choices, task naming, QC, and optional downstream modules. |
| `import.sh` | Raw scanner series classification and acquisition counts; not needed for BIDS import. |
| `tans.yaml` | Optional TANS inputs, hardware, search, and dose settings. |

These shell files are sourced by the launchers, not executed as processing jobs.
Pass the processing config as the **second argument** to `mefmri_pipeline.sh`.
Without that argument, the launcher uses the repository template.

## 3. Set software paths and environments

Edit your copy of `mefmri.sh`:

```bash
EnvironmentScript="/path/to/MyStudy/config/environment.sh"
FS_LICENSE_FILE="/path/to/license.txt"
CHARM_BIN="/path/to/SimNIBS/simnibs_env/bin/charm"  # or "" if charm is on PATH
PIPELINE_PYTHON="/path/to/your/python/environment/bin/python"
```

Set `PIPELINE_PYTHON` near the top, before derived Python settings. It must provide
the internal module dependencies listed in the dependency guide. In
`environment.sh`, configure `FSLDIR`, `FREESURFER_HOME`, Workbench on `PATH`, and
`MSMBINDIR` for MSMSulc. Site environment modules may supply those settings instead.
The example then sources the bundled HCP setup. FreeSurfer needs its own license;
none is distributed here.

The launcher supplies `PFM_MEFMRI_ROOT` before reading the processing config, so
`MEDIR`, priors, module scripts, and other bundled paths keep pointing to the
installation even when the config lives in your study folder. If sourcing the
config directly without the launcher, set `MEDIR="/path/to/pfm-mefmri"` at the top
**before** any `$MEDIR`-derived settings.

`TEDANA_ENV="mefmri_env"`, `WARPKIT_ENV="warpkit_env"`, and `AROMA_ENV="aroma"`
are example environment names, not installed environments. Create or rename them
to match your setup. `*_ACTIVATE_MODE="direct"` uses tools on the current `PATH`;
`conda_activate` and `conda_run` require Conda and the named environment.
`PIPELINE_PYTHON` does not automatically select the tedana, warpkit, or TANS runtime.

For crawling-seed movies, install `wb_surfer2` and FFmpeg. Either put `wb_surfer2`
on `PATH`, set `CRAWLING_SEED_FC_WB_SURFER2` to its executable, or set
`CRAWLING_SEED_FC_WB_SURFER2_CONDA_ENV`. Set `CRAWLING_SEED_FC_ENABLE=0` to skip
that optional movie while keeping other QC.

## 4. Choose processing settings

Review these in your study's `mefmri.sh`:

| Setting | What to choose |
| --- | --- |
| `FUNC_DIRNAME`, `FUNC_FILE_PREFIX` | Match imported task folders and filenames; defaults are `rest` and `Rest`. |
| `PROCESSING_MODE` | `auto`, `multi_echo`, or `single_echo`; auto uses single-echo processing for fewer than three echoes. |
| `DISTORTION_CORRECTION_MODE` | `medic` requires echo-wise phase and warpkit; `topup` requires opposing-PE fieldmaps. Select for your actual acquisition. |
| `MULTI_ECHO_DENOISE_METHOD` | `meica`, `acompcor`, or `aroma`; echo combination still uses tedana for multi-echo branches. |
| `SINGLE_ECHO_DENOISE_METHOD` | `acompcor` or `aroma`. |
| `MEPCA` | Named tedana PCA selection (`kundu` by default), or an intentional component count appropriate for run length. |
| `START_FROM_MODULE`, `STOP_AFTER_MODULE` | Start at `validate`; empty stop runs all enabled stages. |
| `THREADS_*`, `*_PARALLEL_JOBS` | Match available CPU and memory. |
| `CONCAT_FD_THRESHOLD` | Review the censoring threshold for your study. |

The template enables standard QC, concatenation, and NSI. PFM, TANS, and the
experimental multi-network analysis start disabled. To add PFM:

```bash
PFM_ENABLE=1
PFM_STRATEGY="ridgepfm"
PFM_RF_FC_WEIGHT=1.0
PFM_RF_SPATIAL_WEIGHT=0.05
```

Infomap remains available explicitly and requires its separate executable.
For direct-B0 correction, set the fieldmap echo-time difference and filename
patterns for your sequence; the advanced example values are not acquisition
metadata. No-fieldmap processing must also be selected explicitly.

## 5. Configure import and prepare inputs

**Raw DICOM:** edit `import.sh` before running the importer. Set the protocol name,
expected runs, echoes, SBrefs and AP/PA fieldmaps, and the Python regexes matching
your scanner's `SeriesDescription`. The initial `(?!)` regex matches nothing;
retain it only for absent series types. Never replace it with an empty regex,
which matches everything. Use `IMPORT_IGNORE_REGEXES` only for intentional exclusions.

Counts are exact expectations: zero means none expected. Set the expected task
volume count if known; zero disables that volume-count check. Optional minimum
file sizes also use zero to disable checks. The initial anatomy/SBref/fieldmap
volume expectations are one and must be changed for other layouts.
After editing, set `IMPORT_CONFIG_READY=1`; the untouched template stops with
an explanation instead of interpreting your acquisition using another study's rules.

```bash
bash /path/to/pfm-mefmri/bin/mefmri_import.sh \
  /path/to/raw_export /path/to/MyStudy/SUB001 \
  /path/to/MyStudy/config/import.sh --session 1 \
  --import-work-subdir preview-session-1 --dry-run
```

Inspect the import report under `SUB001/import/preview-session-1/`. The dry run
performs conversion and writes working files/reports, but avoids final organized
subject-file staging. Use a fresh working subfolder for the real import:

```bash
bash /path/to/pfm-mefmri/bin/mefmri_import.sh \
  /path/to/raw_export /path/to/MyStudy/SUB001 \
  /path/to/MyStudy/config/import.sh --session 1 --import-work-subdir session-1
```

Use a new working subfolder for each additional preview/import; the importer
refuses to overwrite an existing conversion folder. NORDIC is optional: add `--nordic` only after separately installing its
code and MATLAB and setting `NORDIC_CODE_DIR` or `--nordic-code-dir`.
`mefmri_import_raw_config.sh` is retained as a compatibility template; new studies
should use `mefmri_import_config.sh`.

**BIDS NIfTI:** the BIDS importer uses command-line options rather than the raw
scanner config:

```bash
bash /path/to/pfm-mefmri/bin/mefmri_import.sh \
  /path/to/bids SUB001 /path/to/MyStudy/SUB001 --task rest --mode symlink
```

Use the actual BIDS subject label in place of `SUB001`. See
[input layout and import details](setup-and-inputs.md#importing-inputs).

## 6. Validate one subject, then process

Check shell syntax first:

```bash
bash -n /path/to/MyStudy/config/environment.sh
bash -n /path/to/MyStudy/config/mefmri.sh
```

In `mefmri.sh`, set `START_FROM_MODULE="validate"` and
`STOP_AFTER_MODULE="validate"`, then run:

```bash
bash /path/to/pfm-mefmri/bin/mefmri_pipeline.sh \
  /path/to/MyStudy/SUB001 /path/to/MyStudy/config/mefmri.sh
```

Review the validation output, resolve errors, then set `STOP_AFTER_MODULE=""`
and rerun the same command for the enabled processing chain. Validation checks
inputs; it does not prove that every later module or external tool will succeed.
The subject argument is the subject directory, not the study root.

Keep `RUN_CONFIG_SNAPSHOT=1` to record effective run settings in the subject's QC
outputs. For repeated visits and separate analysis epochs, see
[individual and longitudinal studies](individual-longitudinal-studies.md).

## 7. Add TANS when needed

Follow the [TANS integration guide](tans-integration.md). Copy the YAML template,
replace all software/coil paths and calibration IDs, choose a target, and only
then set `TANS_ENABLE=1`. The shell fragment is appended to the complete processing
config; it is not a replacement for that config. Blank hardware overrides use
the YAML settings. The template makes no automatic 180-degree coil-equivalence
assumption; configure orientation handling for the selected coil.

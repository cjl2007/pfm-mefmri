# Dependencies and installation

Start with a Linux environment and the [configuration guide](configuration.md).
The repository includes pipeline scripts, an HCP Pipelines snapshot, templates,
and a pinned pfm-nsi Git submodule. It does **not** install FSL, FreeSurfer,
Workbench, SimNIBS, Conda, or your Python environments.

The [software inventory](../SOFTWARE_DEPENDENCIES.txt) records the development
stack and module-level commands. These are reference versions, not a guarantee
that every version combination or optional branch has been validated.

## Core software

| Tool | Used for | Installation/source |
| --- | --- | --- |
| FSL | Functional preprocessing, registration, TOPUP, image operations | [FSL installation](https://fsl.fmrib.ox.ac.uk/fsl/docs/install/index.html) |
| FreeSurfer and its license | Anatomy, surfaces, registration | [FreeSurfer installation](https://freesurfer.net/fswiki/DownloadAndInstall) |
| Connectome Workbench | CIFTIs, surface mapping, network outputs | [Workbench downloads](https://www.humanconnectome.org/software/get-connectome-workbench) |
| HCP-compatible MSM | Default `HCP_REGNAME="MSMSulc"` surface registration | Use a build supporting `--regoption`; set `MSMBINDIR`. `HCP_REGNAME="FS"` skips MSMSulc. |
| SimNIBS CHARM | Default CHARM anatomical masking; also used separately by TANS | [SimNIBS installation](https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html) |
| Bash, GNU tools, GNU Parallel, Python | Pipeline orchestration and internal modules | Install through your OS/site package manager; see the software inventory. |

Initialize the NSI dependency with `git submodule update --init --recursive`.
The wrapper uses `lib/pfm-nsi`; its Python dependencies still need installation.

## Python and multi-echo processing

Use a dedicated Python environment with `numpy`, `scipy`, `nibabel`, `pandas`,
`matplotlib`, and `scikit-learn`. Set `PIPELINE_PYTHON` to its interpreter.
`imageio` is needed for optional GIF outputs and `openpyxl` for Excel correction
tables. `Pillow` is used by the optional release-documentation checker.

Multi-echo processing also needs [tedana](https://tedana.readthedocs.io/en/stable/installation.html),
including the `tedana`, `t2smap`, and `ica_reclassify` commands. The recorded
pipeline stack uses Python 3.10 and tedana 26.0.0. For example:

```bash
conda create -n mefmri_env python=3.10 pip
conda run -n mefmri_env python -m pip install \
  'tedana==26.0.0' numpy scipy nibabel pandas matplotlib scikit-learn imageio openpyxl
conda run -n mefmri_env tedana --help
conda run -n mefmri_env ica_reclassify --help
```

This is an environment starting point, not a fully pinned reproduction of all
recorded dependencies. Follow upstream installation instructions when adapting
versions for your platform. The config's example Conda environment names must
match environments you actually create.

## Dependencies selected by processing options

| Option | Additional dependency/setup |
| --- | --- |
| Raw DICOM import | [dcm2niix](https://github.com/rordenlab/dcm2niix); edit the raw-import config for scanner naming and counts. |
| MEDIC (`DISTORTION_CORRECTION_MODE="medic"`) | [warpkit](https://github.com/vanandrew/warpkit), its `wk-*` commands, and echo-wise phase data; use a separate environment such as `warpkit_env`. |
| ICA-AROMA | [ICA-AROMA](https://github.com/maartenmennes/ICA-AROMA); set `AROMA_ENV` or `AROMA_BIN_OVERRIDE`. |
| N4 bias correction (`APPLY_N4_BIAS=1`) | ANTs `N4BiasFieldCorrection`. |
| Crawling-seed FC movie | [wb_surfer2](https://github.com/vanandrew/wbsurfer2), Workbench, and FFmpeg (`ffmpeg`, `ffprobe`); configure its executable/environment or set `CRAWLING_SEED_FC_ENABLE=0`. |
| MEICA QC movies | FFmpeg; optional GIF output also needs `imageio`. |
| RidgePFM | Core Python/Workbench tools and bundled priors; enable with `PFM_ENABLE=1`. |
| Infomap | Separate `infomap` executable; set `PFM_INFOMAP_BINARY` if not on `PATH`. |
| TANS | Separate [pfm-tans installation](https://github.com/cjl2007/pfm-tans), SimNIBS Python runtime, Workbench, coil models and matched calibrations; see [integration](tans-integration.md). |
| NORDIC at import | Separately installed NORDIC code and MATLAB; set `NORDIC_CODE_DIR`. NORDIC source is not bundled. |

MATLAB is not required for the Python PFM/TANS workflow. It is required if you
choose the optional NORDIC import branch. Container recipes are development
starting points and do not imply that every optional dependency is installed.

## Check your environment

After configuring/loading your site environment, use the command checks in the
[software inventory](../SOFTWARE_DEPENDENCIES.txt). At minimum, verify the selected
branches' executables and your Python imports:

```bash
command -v fslmaths recon-all wb_command parallel charm
/path/to/your/python/environment/bin/python -c \
  'import numpy, scipy, nibabel, pandas, matplotlib, sklearn'
```

Also check tedana, warpkit, or AROMA within their own configured environments.
Then run [one subject through validation](configuration.md#6-validate-one-subject-then-process)
before starting the full pipeline. No package installer in this guide launches
participant processing.

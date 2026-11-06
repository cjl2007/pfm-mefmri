# pfm-mefmri

Processing and denoising for individual resting-state and task fMRI, from a single
visit to densely sampled longitudinal studies. The workflow supports multi-echo
and single-echo data, multiple distortion-correction strategies, careful
coregistration to anatomy, and consistent within-participant alignment over time.
Quality control connects anatomical alignment, motion and physiology, expert-calibrated FC quality scores (Lynch et al., 2026), and
crawling-seed FC movies. Optional modules add precision functional mapping and
individualized TMS targeting with [TANS](https://github.com/cjl2007/pfm-tans).

![The dataset journey through pfm-mefmri: HCP anatomical pipelines, functional correction, longitudinal alignment, denoising, CIFTI outputs, QC, PFM, and individualized TMS targeting](docs/assets/pipeline-journey.svg)

**[Open the interactive workflow and network-map viewers](https://cjl2007.github.io/pfm-mefmri/)**

## Explore the workflow

| Section | What you will find |
| --- | --- |
| [Data organization](docs/data-organization.md) | Expected inputs, validate checks, and example rest/Motor/WorkingMemory outputs. |
| [Processing and setup](docs/pipeline-overview.md) | Distortion correction, denoising, anatomy, longitudinal alignment, and CIFTI outputs. |
| [Within-subject longitudinal alignment](docs/longitudinal-alignment.md) | Full-sequence comparison followed by focused examples of alignment within the same participant across scans and sessions. |
| [Individual network maps](#precision-functional-mapping-pfm) | RidgePFM WTA and softmax views across network taxonomies; conventional Infomap is also supported. |
| [RidgePFM parameters](docs/ridgepfm-parameters.md) | Interactive ME01 maps followed by reliability–individuality trade-offs and practical guidance. |
| [Quality-control gallery](docs/output-gallery.md) | Motion and physiology, expert-calibrated FC quality scores, crawling-seed FC movies, alignment, and signal QC. |
| [Individualized TMS](docs/tans/README.md) | Basic workflow, coil-search movie, and optional TANS module setup. |

## Within-subject longitudinal alignment

[Review alignment across scans and sessions](docs/longitudinal-alignment.md).
Both approaches align well to anatomy; within-subject refinement adds subtle
fine-tuning to dampen occasional placement shifts over time. This attention to
within-person consistency is a central design priority of the pipeline.

## Precision Functional Mapping (PFM)

**RidgePFM** combines individual functional connectivity with canonical network
templates and spatial priors whose [influence you can weight and explore](docs/ridgepfm-parameters.md).
It is a work in progress. A conventional **Infomap** workflow is also supported
as an option.

[![Individual network assignments; open the viewer to select a network taxonomy and any softmax network](docs/assets/ridgepfm/taxonomy-explorer-poster.png)](https://raw.githubusercontent.com/cjl2007/pfm-mefmri/main/docs/ridgepfm-taxonomy-explorer.zip)

**[Open the network-map selector](https://raw.githubusercontent.com/cjl2007/pfm-mefmri/main/docs/ridgepfm-taxonomy-explorer.zip)**:
choose **Lynch2024**, **Du2024**, **Hermosillo2024**, or **Kong2019**, then WTA
labels or any network’s full-data softmax map. Check **Areal parcellation** to
view parcels derived from that taxonomy’s network map and individual FC.
Parcel fills show the majority network identity or mean softmax support;
borders follow the selected taxonomy’s parcels. With the option off, network
borders remain visible on softmax maps. Low–High colors show relative network
support. Download and extract the ZIP, then open the HTML in a browser.

[Explore λ and spatial-prior effects, then view the trade-off charts](docs/ridgepfm-parameters.md) ·
[Mapping outputs](docs/pfm-outputs.md) · [Sampling across visits](docs/individual-longitudinal-studies.md).

## Data organization

[Expected inputs and example processed outputs](docs/data-organization.md) for a
participant with resting-state, Motor, and WorkingMemory scans across sessions
and runs, including the naming and metadata checked by `validate`.

## Get started

1. Install the [dependencies](SOFTWARE_DEPENDENCIES.txt).
2. Adapt the [study configuration](config/mefmri_wrapper_config.sh).
3. [Import your data](docs/setup-and-inputs.md#importing-inputs), then run:

```bash
bash /path/to/pfm-mefmri/bin/mefmri_pipeline.sh \
  /data/MyStudy/SUB001 /data/MyStudy/config/mefmri.sh
```

[Quick start](docs/quick-start.md) · [Workflow stages](docs/pipeline-overview.md) ·
[Output locations](docs/output-map.md) · [All guides](docs/README.md).

Research software. Review participant data quality and modeled outputs before
use; participant-derived figures require permission before public distribution.

[Configure your project](docs/configuration.md): copy the generic templates outside the repository, set software paths and acquisition rules, and choose your processing tasks.

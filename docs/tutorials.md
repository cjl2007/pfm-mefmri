# Tutorials

These tutorials are written as practical workflows. They assume the input data
are already organized in the expected subject directory layout. First complete
[configuration setup](configuration.md); edit the settings shown below in that
study config and pass its path to each pipeline command.

## Tutorial 1: Run One Resting-State Subject

Goal: produce preprocessing, QC, concatenated CIFTI, NSI, and PFM outputs for
one subject.

Config choices:

```bash
FUNC_DIRNAME="rest"
FUNC_FILE_PREFIX="Rest"
START_FROM_MODULE="validate"
STOP_AFTER_MODULE=""
PROCESSING_MODE="auto"
DISTORTION_CORRECTION_MODE="medic"
MULTI_ECHO_DENOISE_METHOD="meica"
CONCAT_ENABLE=1
NSI_ENABLE=1
PFM_ENABLE=1
PFM_STRATEGY="ridgepfm"
```

Run:

```bash
bash bin/mefmri_pipeline.sh /path/to/study/ME06 /path/to/MyStudy/config/mefmri.sh
```

Review in this order:

1. `func/rest/QC/Validation/`
2. `func/rest/QC/Coreg/functional_coregistration_qc.html`
3. `func/rest/QC/Fieldmaps/distortion_correction_qc.html`
4. `func/rest/QC/MotionGrayPlots/motion_grayplots_qc.html`
5. `func/rest/QC/NSI/`
6. `func/rest/PFM/RidgePFM/RidgePFM_VTX.dlabel.nii`

Figure slots:

```text
docs/assets/tutorial-one-subject-coreg.png
docs/assets/tutorial-one-subject-motion.png
docs/assets/tutorial-one-subject-ridgefusion.png
```

## Tutorial 2: Rerun Only QC Reports

Goal: regenerate reports after preprocessing outputs already exist.

Coreg/SDC QC:

```bash
COREG_SDC_QC_ENABLE=1
START_FROM_MODULE="coreg"
STOP_AFTER_MODULE="coreg"
```

Motion grayplots QC:

```bash
MOTION_GRAYPLOTS_QC_ENABLE=1
START_FROM_MODULE="mgtr"
STOP_AFTER_MODULE="mgtr"
```

The exact rerun window depends on which outputs already exist. The QC scripts
read existing preprocessing outputs and write HTML/reportlet files under
`func/<FUNC_DIRNAME>/QC/`.

## Tutorial 3: Check Whether A Dataset Is Usable For PFM

Goal: run NSI usability on the FD-censored concatenated CIFTI.

Config:

```bash
NSI_ENABLE=1
NSI_USE_EXTERNAL_CLI=1
NSI_USABILITY_MODEL=1
NSI_RELIABILITY_MODEL=0
START_FROM_MODULE="nsi"
STOP_AFTER_MODULE="nsi"
```

Run:

```bash
bash bin/mefmri_pipeline.sh /path/to/study/ME06 /path/to/MyStudy/config/mefmri.sh
```

Inspect:

```text
func/rest/QC/NSI/*_run_config.txt
func/rest/QC/NSI/*_nsi.npz
func/rest/QC/NSI/*_hist_nsi.png
func/rest/QC/NSI/*_hist_nsi_by_network.png
func/rest/QC/NSI/*_hist_nsi_by_structure.png
```

Interpretation:

- Broadly higher NSI supports better expected PFM usability.
- Network/structure histograms reveal whether quality is system-specific.
- The run config confirms the input CIFTI, censoring threshold, duration,
  structures, memory controls, and exact CLI call.

Figure slots:

```text
docs/assets/tutorial-nsi-usability-hist.png
docs/assets/tutorial-nsi-network-hist.png
```

## Tutorial 4: Estimate Reliability From A Scout Scan

Goal: ask how reliable a longer map may be from a shorter observed sample.

Config:

```bash
NSI_ENABLE=1
NSI_USE_EXTERNAL_CLI=1
NSI_USABILITY_MODEL=0
NSI_RELIABILITY_MODEL=1
NSI_RELIABILITY_NSI_T=10
NSI_RELIABILITY_QUERY_T=60
NSI_EXTERNAL_THRESHOLDS="0.6,0.7,0.8"
START_FROM_MODULE="nsi"
STOP_AFTER_MODULE="nsi"
```

Inspect:

```text
func/rest/QC/NSI/*_reliability.json
func/rest/QC/NSI/*_run_config.txt
```

Use this when deciding whether a participant needs more scan time before PFM.

## Tutorial 5: Review RidgePFM PFM

Goal: inspect the default fast PFM path.

Config:

```bash
PFM_ENABLE=1
PFM_STRATEGY="ridgepfm"
PFM_HOMOGENEITY_TEST_ENABLE=1
PFM_AREAL_ENABLE=0
START_FROM_MODULE="pfm"
STOP_AFTER_MODULE="pfm"
```

Inspect:

```text
func/rest/PFM/RidgePFM/RidgePFM_VTX.dlabel.nii
func/rest/PFM/RidgePFM/RidgePFM_VTX_R2.dtseries.nii
func/rest/PFM/RidgePFM/RidgePFM_VTX_ProbMaps.dtseries.nii
func/rest/PFM/RidgePFM/HomogeneityTest/RidgePFMHomogeneity_summary.png
```

In Workbench:

- Load the subject midthickness or inflated surfaces.
- Load `RidgePFM_VTX.dlabel.nii` as the primary label map.
- Load probability/R2 maps as supporting layers.
- Compare the map against QC reports before interpreting unusual network
  boundaries.

Figure slots:

```text
docs/assets/tutorial-ridgefusion-workbench.png
docs/assets/tutorial-ridgefusion-homogeneity.png
```

## Tutorial 6: Review Infomap And Apply Manual Corrections

Goal: run graph-density Infomap, review algorithmic labels, then apply manual
corrections without rerunning Infomap.

Config for Infomap:

```bash
PFM_ENABLE=1
PFM_STRATEGY="infomap"
PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1
PFM_HOMOGENEITY_TEST_ENABLE=1
START_FROM_MODULE="pfm"
STOP_AFTER_MODULE="pfm"
```

Inspect:

```text
func/rest/PFM/Infomap/GraphDensity_<value>/
func/rest/PFM/Infomap/InfomapNetworkLabels_ModeConsensus.dlabel.nii
func/rest/PFM/Infomap/InfomapNetworkLabels_ProbabilityConsensus.dscalar.nii
func/rest/PFM/Infomap/HomogeneityTest/InfomapDensityHomogeneity_summary.png
```

Edit the per-density manual correction tables:

```text
GraphDensity_<value>/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv
```

Then rerun only the update stage:

```bash
PFM_INFOMAP_UPDATE_ENABLE=1
START_FROM_MODULE="pfm_update"
STOP_AFTER_MODULE="pfm_update"
```

Run:

```bash
bash bin/mefmri_pipeline.sh /path/to/study/ME06 /path/to/MyStudy/config/mefmri.sh
```

Inspect the adjusted consensus outputs:

```text
func/rest/PFM/Infomap/InfomapNetworkLabels_ManualAdjusted_ModeConsensus.dlabel.nii
func/rest/PFM/Infomap/InfomapNetworkLabels_ManualAdjusted_ProbabilityConsensus.dscalar.nii
```

Figure slots:

```text
docs/assets/tutorial-infomap-density-review.png
docs/assets/tutorial-infomap-manual-corrections.png
docs/assets/tutorial-infomap-adjusted-consensus.png
```

## Tutorial 7: Build The Crawling Seed FC Movie

Goal: create a fast visual scan of seed-based FC structure from the final
concatenated CIFTI.

Config:

```bash
CRAWLING_SEED_FC_ENABLE=1
START_FROM_MODULE="concat"
STOP_AFTER_MODULE="concat"
```

Run:

```bash
bash bin/mefmri_pipeline.sh /path/to/study/ME06 /path/to/MyStudy/config/mefmri.sh
```

Inspect:

```text
func/rest/QC/CrawlingSeedFC/sub-ME06.mp4
```

This module can be slow and memory-intensive because it creates a dense
correlation CIFTI before rendering. For routine batches, disable it unless the
movie is needed:

```bash
CRAWLING_SEED_FC_ENABLE=0
```

Figure/movie slots:

```text
docs/assets/tutorial-crawling-seed-fc-still.png
docs/assets/tutorial-crawling-seed-fc-demo.mp4
```

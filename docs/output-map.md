# Output Map

This page lists the output locations most users inspect after a run. The exact
files depend on config settings, data type, and enabled modules.

See the [illustrated output gallery](output-gallery.md) for example reportlets
and interpretation, including the optional TANS handoff.

## TANS (optional)

```text
<SubjectDir>/tans/
  reports/index.html
  targets/<target>/reports/
  shared/export_<subject>/index.html
  shared/export_<subject>/<target>/dose-calibration/DoseCalibrationReport.html
```

`TANS_OUTPUT_DIRNAME` changes the subject-relative root. Preserve the complete
operator package; save downloaded session records outside it. See
[TANS integration](tans-integration.md).

## Functional QC Root

Controlled by:

```bash
QC_OUTPUT_DIRNAME="QC"
```

Default location:

```text
<SubjectDir>/func/<FUNC_DIRNAME>/QC/
```

## Run Metadata

```text
func/<FUNC_DIRNAME>/QC/RunMetadata/pipeline_run_<timestamp>.txt
```

Written when:

```bash
RUN_CONFIG_SNAPSHOT=1
```

## Coregistration And Fieldmap QC

```text
func/<FUNC_DIRNAME>/QC/Coreg/
  functional_coregistration_qc.html
  coreg_stability_metrics.tsv
  manifest.json
  reportlets/

func/<FUNC_DIRNAME>/QC/Fieldmaps/
  distortion_correction_qc.html
  run_metrics.tsv
  manifest.json
  reportlets/
```

## Motion Grayplots

```text
func/<FUNC_DIRNAME>/QC/MotionGrayPlots/
  motion_grayplots_qc.html
  motion_grayplot_metrics.tsv
  manifest.json
  reportlets/
```

## Concatenated CIFTIs

```text
func/<FUNC_DIRNAME>/ConcatenatedCiftis/
  <FuncPrefix>_<InputTag>_Concatenated.dtseries.nii
  <FuncPrefix>_<InputTag>_Concatenated+FDlt<FD>.dtseries.nii
  ScanInfo.json
```

The FD-censored file is the default input for NSI, crawling seed FC movies, and
PFM.

## NSI

```text
func/<FUNC_DIRNAME>/QC/NSI/
  <prefix>_run_config.txt
  <prefix>_nsi.npz
  <prefix>_hist_nsi.png
  <prefix>_hist_nsi_by_network.png
  <prefix>_network_hist_summary.csv
  <prefix>_hist_nsi_by_structure.png
  <prefix>_structure_hist_summary.csv
  <prefix>_reliability.json
```

The default prefix is derived from usability/reliability settings and total
concatenated duration.

## Crawling Seed FC

```text
func/<FUNC_DIRNAME>/QC/CrawlingSeedFC/
  sub-<Subject>.mp4
```

Optional/intermediate files may include rewritten scenes, dtseries links, source
movies, and dense correlation files depending on `CRAWLING_SEED_FC_*` settings.

## PFM RidgePFM

```text
func/<FUNC_DIRNAME>/PFM/RidgePFM/
  RidgePFM_VTX.dlabel.nii
  RidgePFM_VTX.L.border
  RidgePFM_VTX.R.border
  RidgePFM_VTX_TemplateFit_R2.dtseries.nii
  RidgePFM_VTX_FusedScoreSoftmax.dtseries.nii
  RidgePFM_VTX+ArealParcellation.dlabel.nii
  RidgePFM_VTX+ArealParcellation.L.border
  RidgePFM_VTX+ArealParcellation.R.border
  RidgePFM_VTX+ArealParcellation_FusedScoreSoftmax.dtseries.nii
  RepeatAssignment/
  RidgePFM_Vertexwise_Review.scene
  RidgePFM_ArealParcellation_Review.scene
  prep/
  HomogeneityTest/
```

## PFM Infomap

```text
func/<FUNC_DIRNAME>/PFM/Infomap/
  Bipartite_PhysicalCommunities.dtseries.nii
  InfomapNetworkLabels_ModeConsensus.dlabel.nii
  InfomapNetworkLabels_ProbabilityConsensus.dscalar.nii
  InfomapNetworkLabels_ManualCorrections.csv
  InfomapNetworkLabels_ManualCorrections.xlsx
  GraphDensity_<value>/
  HomogeneityTest/
```

## Areal Parcellation

```text
func/<FUNC_DIRNAME>/PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.dlabel.nii
func/<FUNC_DIRNAME>/PFM/Infomap/InfomapNetworkLabels+ArealParcellation.dlabel.nii
```

The exact prefix depends on strategy and manual-label settings.

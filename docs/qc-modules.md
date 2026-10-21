# QC Modules

The [output gallery](output-gallery.md) includes actual coregistration,
distortion-correction, motion, NSI, and seed-FC examples. Figure slots below are
suggested report families; only the explicitly linked gallery assets are bundled.

The pipeline includes multiple QC layers. They are intended to answer different
questions: whether inputs were understood, whether anatomy and field maps line
up, whether motion remains coupled to signal, whether ME-ICA component decisions
look plausible, whether the data support PFM, and whether final network maps are
stable enough to interpret.

The functional QC root is controlled by:

```bash
QC_OUTPUT_DIRNAME="QC"
```

Most QC outputs therefore live under:

```text
func/<FUNC_DIRNAME>/QC/
```

## Input Validation

Enable with:

```bash
QC_VALIDATION_ENABLE=1
```

Use this first when onboarding a new dataset. It checks that the expected task,
session, run, echo, SBRef, JSON, and fieldmap inputs are visible before the
heavier stages begin.

Expected location:

```text
func/<FUNC_DIRNAME>/QC/Validation/
```

Figure slot:

```text
docs/assets/qc-validation-summary.png
```

## Coregistration And SDC QC

Enable with:

```bash
COREG_SDC_QC_ENABLE=1
COREG_SDC_QC_STRICT=0
```

Set strict mode to `1` only when a QC generation failure should abort the
pipeline.

Primary reports:

```text
func/<FUNC_DIRNAME>/QC/Coreg/functional_coregistration_qc.html
func/<FUNC_DIRNAME>/QC/Fieldmaps/distortion_correction_qc.html
```

Common reportlets include:

```text
Coreg/reportlets/subj_desc-coreg_T1wAnchor.png
Coreg/reportlets/subj_desc-coreg_unwarpedSBRefBBR_fade.gif
Coreg/reportlets/subj_desc-coreg_runSimilarityMatrix.png
Coreg/reportlets/subj_desc-coreg_anatSimilarityOverTime.png
Fieldmaps/reportlets/subj_desc-sdc_<Sx_Ry>_flicker.gif
Fieldmaps/run_metrics.tsv
Coreg/coreg_stability_metrics.tsv
```

Use these reports to inspect:

- Whether SBRefs align to the T1w anatomy.
- Whether fieldmap or MEDIC correction improves the distorted reference.
- Whether run-to-run registration is stable across sessions.
- Whether the same anatomical region is sampled consistently over time.

Figure slots:

```text
docs/assets/qc-coreg-bbr.png
docs/assets/qc-sdc-flicker.gif
docs/assets/qc-coreg-stability.png
```

## Motion Grayplots And QC-FC Proxy

Enable with:

```bash
MOTION_GRAYPLOTS_QC_ENABLE=1
MOTION_GRAYPLOTS_QC_STRICT=0
MOTION_GRAYPLOTS_QC_FD_THRESHOLD=0.3
```

Primary report:

```text
func/<FUNC_DIRNAME>/QC/MotionGrayPlots/motion_grayplots_qc.html
```

Important files:

```text
MotionGrayPlots/reportlets/subj_desc-motionGrayplots_stacked.png
MotionGrayPlots/reportlets/subj_desc-motionCouplingSummary.png
MotionGrayPlots/reportlets/subj_desc-motionGrayplot_S<session>_R<run>.png
MotionGrayPlots/motion_grayplot_metrics.tsv
MotionGrayPlots/manifest.json
```

The report stacks multiple processing stages and overlays motion, FD,
DVARS, FD censor markers, and run boundaries. When Siemens RESP/PULS logs and
their acquisition-info log are present under
`physio/unprocessed/<func>/session_<N>/run_<N>/`, it also aligns respiratory and
pulse traces to scanner acquisition ticks and displays them above FD. Missing
physiology remains optional and does not prevent report generation. The coupling
summary reports motion/signal relationships before and after censoring.

FD uses the same historical calculation as [CIFTI concatenation](pfm-outputs.md):
motion parameters are bandstop-filtered only for TR < 2 seconds, from 0.2 Hz
to Nyquist minus 0.019 Hz. Frames at or above the threshold are censored.
Set `MOTION_GRAYPLOTS_QC_FD_THRESHOLD` equal to `CONCAT_FD_THRESHOLD` for matching
retained-frame counts. Invalid or too-short motion inputs fail explicitly
instead of silently switching to unfiltered FD.

Review priorities:

- Look for run-specific signal disruptions.
- Compare pre-denoising and post-denoising stages.
- Check whether high-motion frames dominate DVARS or grayordinate structure.
- Inspect `median_abs_grayordinate_fd_r` and `p95_abs_grayordinate_fd_r` in the
  metrics TSV for residual motion coupling.

Figure slots:

```text
docs/assets/qc-motion-grayplots-stacked.png
docs/assets/qc-motion-coupling-summary.png
```

## ME-ICA Reclassification QC

The ME-ICA branch supports component reclassification and report generation:

```bash
MEICA_RECLASSIFY_ENABLE=1
MEICA_CLASSIFIER_MODE="nsi"       # nsi|legacy_template_rho|none
MEICA_RECLASSIFY_IMPLEMENTATION="standard"
MEICA_RECLASS_NO_REPORTS=0
MEICA_QC_CIFTI_ENABLE=1
```

The production default is the standard NSI reclassifier.

Experimental spatial/tissue QC outputs are kept as developer/testing options
rather than normal user-facing denoising settings:

```bash
MEICA_RECLASSIFY_SPATIAL_QC_MOVIES=1
MEICA_RECLASSIFY_SPATIAL_QC_MOVIE_MP4=1
MEICA_RECLASSIFY_SPATIAL_QC_MOVIE_GIF=0
```

Common files include:

```text
ComponentDecisions.tsv
ClassificationSummary.json
ClassificationSummary.txt
KappaRhoLog_vs_NSI.png
KappaRhoLog_vs_NSI_spatial_traffic.png
Predicted_vs_Observed_LogKappaRho.png
qc_media_summary.json
index.html
movies/AcceptedFinal.mp4
movies/RejectedFinal.mp4
movies/Killed.mp4
movies/Rescued.mp4
```

Use this QC to inspect which components were accepted, rejected, killed,
rescued, or marked for review, and whether classification decisions are
consistent with kappa/rho/NSI evidence.

Figure slots:

```text
docs/assets/qc-meica-kappa-rho-nsi.png
docs/assets/qc-meica-spatial-traffic.png
docs/assets/qc-meica-component-movie-still.png
```

## NSI Usability And Reliability

Enable with:

```bash
NSI_ENABLE=1
NSI_USE_EXTERNAL_CLI=1
NSI_USABILITY_MODEL=1
NSI_RELIABILITY_MODEL=0
```

NSI uses the concatenated, FD-censored CIFTI by default and writes to:

```text
func/<FUNC_DIRNAME>/QC/NSI/
```

Common outputs:

```text
pfm_nsi_usability_full<minutes>m_run_config.txt
pfm_nsi_usability_full<minutes>m_nsi.npz
pfm_nsi_usability_full<minutes>m_hist_nsi.png
pfm_nsi_usability_full<minutes>m_hist_nsi_by_network.png
pfm_nsi_usability_full<minutes>m_hist_nsi_by_structure.png
pfm_nsi_usability_full<minutes>m_network_hist_summary.csv
pfm_nsi_usability_full<minutes>m_structure_hist_summary.csv
pfm_nsi_reliability_nsiT<minutes>m_reliability.json
```

Interpretation:

- Usability summarizes whether the observed data quality and duration are likely
  adequate for PFM.
- Reliability estimates expected map reliability at a requested duration.
- Network and structure histograms show whether NSI is uneven across systems or
  anatomical structures.
- The run config records the exact CLI arguments and resource controls used.

Figure slots:

```text
docs/assets/qc-nsi-histogram.png
docs/assets/qc-nsi-network-histograms.png
docs/assets/qc-nsi-usability-curve.png
```

## PFM Homogeneity Diagnostics

Enable with:

```bash
PFM_HOMOGENEITY_TEST_ENABLE=1
PFM_HOMOGENEITY_ROTATIONS_CIFTI="$MEDIR/res0urces/Rotated_inds.dtseries.nii"
```

Outputs are written under the active PFM strategy:

```text
func/<FUNC_DIRNAME>/PFM/RidgePFM/HomogeneityTest/
func/<FUNC_DIRNAME>/PFM/Infomap/HomogeneityTest/
```

Common files:

```text
RidgePFMHomogeneity_summary.csv
RidgePFMHomogeneity_summary.png
RidgePFMHomogeneity_nulls.npy
InfomapDensityHomogeneity_summary.csv
InfomapDensityHomogeneity_summary.png
InfomapDensityHomogeneity_nulls.npy
```

For Infomap, this diagnostic is especially useful across graph densities. When
spin rotations are available, the figure can show empirical homogeneity against
spin-null expectations and mark the selected density.

Figure slot:

```text
docs/assets/qc-pfm-homogeneity-summary.png
```

## Crawling Seed FC Movie

Enable with:

```bash
CRAWLING_SEED_FC_ENABLE=1
```

This concat-stage QC sub-step runs after concat and before NSI when enabled. To
rerun concat and the optional movie, use:

```bash
START_FROM_MODULE="concat"
STOP_AFTER_MODULE="concat"
```

Primary output:

```text
func/<FUNC_DIRNAME>/QC/CrawlingSeedFC/sub-<Subject>.mp4
```

The module builds a dense correlation CIFTI and renders a Workbench/wb_surfer2
movie from sampled seed vertices. It is intentionally resource-heavy. The dense
correlation intermediate can be tens of GB and is deleted by default unless:

```bash
CRAWLING_SEED_FC_KEEP_DCONN=1
```

Use this movie to look for broad spatial artifacts, poor surface alignment,
unusual connectivity structure, or session/run patterns that are hard to see in
static plots.

Figure/movie slots:

```text
docs/assets/qc-crawling-seed-fc-still.png
docs/assets/qc-crawling-seed-fc-demo.mp4
```

# Output gallery

The QC, seed-FC, and individual-network images come from one densely sampled
participant with **13 sessions and 26 resting-state runs**. These replace the
earlier single-session layout examples. The TANS illustrations retain a separate
SAL example and do not represent stimulation planning for the QC participant.

This gallery follows the participant from data-quality review to individual
mapping and optional stimulation planning. It contains selected examples, not
automatic pass/fail decisions, treatment recommendations, or a benchmark of
clinical performance.

## Alignment and fieldmap/SDC review

![Animated anatomical anchor and functional-reference alignment across 26 runs](assets/qc-coreg-over-time.gif)

The animation cycles through the T1w anchor and all 26 saved run panels in numeric
session/run order in a fast **8.87-second loop** (220 ms per run panel with
short crossfades). It uses the production fade routine with fewer transition
frames for a compact documentation asset; no registration is rerun.

![One corrected run reference with anatomical outlines in three orthogonal views](assets/qc-coreg-bbr.png)

![Original susceptibility-distortion correction flicker comparison](assets/qc-sdc-flicker.gif)

These are production reportlets from `QC/Coreg/` and `QC/Fieldmaps/`. The
animated panel shows the correction comparison, not an electric-field simulation.
Review the complete HTML reports and metrics, including regions with weak signal
or potential susceptibility-related dropout.

![Anatomical similarity and between-run reference consistency in scan order](assets/qc-coreg-stability.png)

The alignment summary includes **all 26 runs across 13 sessions**. It replots
the saved metric values in numeric session/run order; the original report used
lexicographic folder order. No registration or similarity estimation was rerun.
Repeated-session alignment QC does not measure longitudinal network change or
an intervention effect.

## Signal quality and network signal

![Motion traces, FD, DVARS, and grayplots through denoising and MGTR](assets/qc-motion-grayplots.png)

The motion figure shows **session 1, run 1**, selected by acquisition order for
readability. It uses saved processing-stage CIFTIs, the production QC plotting
code, 1,200 sampled grayordinates, and a generic display title. No denoising is
rerun. Its 0.3 mm FD line reflects the example configuration, not a recommended
universal threshold. The full subject report retains the other runs and the
wide concatenated grayplot.

![Censoring burden and signal-motion coupling across 26 runs in numeric session order](assets/qc-motion-coupling.png)

This unchanged production reportlet summarizes all 26 runs. The lower panels
compare absolute signal-motion correlations through processing stages after FD
censoring; faint dashed traces show the corresponding values before censoring.
Lower coupling alone does not establish that every source of artifact is removed.

![NSI model-estimated PFM usability](assets/qc-nsi-usability.png)

![NSI distributions across example networks](assets/qc-nsi-networks.png)

The curve estimates PFM usability using the selected NSI model. It is not a
direct test–retest measurement, and its operating characteristics should not be
assumed to transfer to every acquisition or population. Network-wise summaries
help reveal heterogeneous quality that a whole-brain average can hide.

## Seed-FC review

![Animated preview of the crawling cortical seed-FC movie](assets/qc-seed-fc.gif)

[Full original seed-FC movie](assets/qc-seed-fc.mp4). The looping GIF samples one
frame per second from the saved movie and plays at four frames per second,
covering the full sweep in approximately 17.5 seconds. The original MP4 retains
all frames. The [still image](assets/qc-seed-fc.png) is decoded at 12 seconds;
no new connectivity maps are computed.
The movie and scene are written to `QC/CrawlingSeedFC/` when enabled.

## Individual network representations

![RidgePFM WTA labels and one continuous network representation](assets/pfm-ridge-fusion.png)

RidgePFM: all-network winner-take-all (WTA) maps on the left; the probabilistic /
weighted map for **DMN-Parietal** (`Default_Parietal`, the first network) on the
right. These documentation views are rendered with Connectome Workbench
from the saved CIFTIs and participant surfaces. The weight output is normalized
model evidence; it is not a calibrated posterior probability of true membership.

![Infomap mode-consensus labels and one network's cross-density frequency](assets/pfm-infomap.png)

Infomap: mode-consensus labels on the left, first network's assignment frequency
on the right. Frequency is computed over densities with a valid assigned label
at that location; unassigned labels are excluded from the denominator. It is
not independent-run or longitudinal reliability. These panels use the same
participant as the QC examples but are not a formal accuracy comparison between
mapping methods. Surface boundaries in the RidgePFM view are its own saved
network borders; no RidgePFM borders are superimposed on the Infomap view.
Each continuous panel uses its stored display palette; color intensity is not
a common quantitative scale for comparing the two methods.

These are saved pooled-data maps, not separately estimated visit maps. The
RidgePFM display uses results generated under the earlier Ridge-Fusion name;
this documentation update does not recompute them with the new spatial weight
of 0.05. Figure provenance and filenames therefore retain the historical label. The
renderings show cortical views; the underlying CIFTIs can also contain
subcortical data. The existing Workbench review scenes and detailed
[PFM output guide](pfm-outputs.md) retain those broader review paths.

## Optional stimulation planning

![Curated TANS optimization example](assets/tans-sal-optimization.png)

![Illustrative technician-facing TANS dose-workpage walkthrough](assets/tans-sal-dose-workpage.png)

These reuse the reviewed SAL examples from pfm-tans. SAL/DLPFC is illustrative:
network identity and cortical search zone are separate choices. The optimization
layout has schematic history scores and an adjacent illustrated candidate rather
than the optimizer's winning placement. The dose-workpage illustration uses a
hypothetical 50% MSO rMT and 55% treatment MSO; the model's selection remains
distinct from the displayed/recorded example intensity.

These curated publication layouts are not guaranteed end-of-run figure layouts.
The actual software generates research reports, search replay outputs when
configured, and matched operator pages. See [TANS integration](tans-integration.md)
and the [pfm-tans gallery](https://github.com/cjl2007/pfm-tans).

## Provenance and sharing

`assets/manifest.json` records each bundled file's fingerprint and rendering
type. Restricted local provenance records retain the source paths. Raw MRI,
participant CIFTIs, head meshes, actual dose/session records, and operator export
packages are not included in this documentation gallery.

All participant-derived figures require maintainer permission review before
public distribution. Generic filenames/titles are not a guarantee of
anonymization. No model output shown here establishes clinical safety, treatment
response, or a universal neuronal activation threshold.

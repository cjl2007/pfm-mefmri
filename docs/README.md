# PFM-ME fMRI Pipeline Documentation

[Current local release update (2026-10-13)](release-2026-10-13.md).

This documentation describes the PFM-ME fMRI processing pipeline as a modular
workflow for multi-echo and single-echo fMRI preprocessing, quality control,
concatenation, Network Similarity Index evaluation, and precision functional
mapping.

The pages are written in a wiki style so they can be used directly in the repo,
copied into a GitHub Wiki, or linked from a shorter top-level README.

## Start Here

- [Quick Start](quick-start.md): minimal setup, expected input layout, and a
  first run.
- [Pipeline Overview](pipeline-overview.md): processing stages, module order,
  and the main configuration switches.
- [QC Modules](qc-modules.md): visual reports, grayplots, NSI, ME-ICA
  reclassification, homogeneity diagnostics, and crawling seed FC movies.
- [PFM Outputs](pfm-outputs.md): RidgePFM, Infomap, areal parcellation, and
  files to inspect in Connectome Workbench.
- [Output Map](output-map.md): where important outputs are written on disk.
- [Tutorials](tutorials.md): practical walkthroughs for preprocessing, QC review,
  NSI, and PFM review.
- [Individual and longitudinal studies](individual-longitudinal-studies.md):
  analysis epochs, pooled references, repeated observations, and intervention provenance.
- [TANS integration](tans/README.md): optional final-module setup, configuration, and input selection.
- [Output gallery](output-gallery.md): illustrated QC, mapping, and stimulation outputs.
- [Local release draft](release-2026-09-21.md): included changes and publication checks.

## What The Pipeline Does

PFM-ME fMRI is designed for dense, longitudinal single-participant processing.
The pipeline keeps run-level provenance while producing subject-level CIFTI and
PFM products that can be inspected across sessions, runs, denoising branches,
QC modules, and network-mapping strategies.

Core capabilities include:

- HCP-style anatomical processing with optional CHARM-derived brain and cortical
  ribbon masks.
- Task-local functional processing, including `FUNC_DIRNAME="All"` for looping
  over multiple task folders.
- Distortion correction through TOPUP, direct B0 field maps, MEDIC/warpkit, or
  explicit no-fieldmap mode.
- Multi-echo processing with tedana/ME-ICA, OCME plus aCompCor, or OCME plus
  ICA-AROMA.
- Single-echo processing with aCompCor or ICA-AROMA.
- Optional mean gray-timecourse regression (MGTR), volume-to-surface mapping,
  CIFTI concatenation, and FD censoring.
- Downstream NSI usability/reliability scoring.
- Precision functional mapping through RidgePFM or Infomap.
- Optional Infomap community-to-network labeling, manual label updates,
  homogeneity diagnostics, and areal parcellation.
- Rich QC outputs intended for fast review before downstream interpretation.

## Figure Assets

The [output gallery](output-gallery.md) now includes representative QC reports,
seed FC, individual maps, continuous representations, and TANS outputs.
[Asset notes](assets/README.md) distinguish unchanged reportlets, display-only
re-renderings, and curated TANS illustrations. Participant-derived assets need
sharing-permission review before public distribution.

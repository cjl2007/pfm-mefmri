# Individual and longitudinal studies

The workflow organizes processing around one participant and retains session/run
provenance. That organization supports dense individual mapping and repeated
observations, but the scientific comparison still needs a study design.

## Decide what each map represents

| Aim | Analysis unit | Important distinction |
| --- | --- | --- |
| Stable individual reference | Pooled, QC-approved runs within a defined sampling period | More data can support a reference map; pooling removes temporal resolution. |
| Change across visits or intervention epochs | Separate, explicitly defined run sets and output roots | Use consistent methods and evaluate reliability before interpreting differences. |
| Stimulation planning | A reviewed map and anatomy selected for a particular planning time point | Record which visit/epoch and software/configuration produced the target. |

Do not repeatedly overwrite one PFM folder and treat it as a longitudinal record.
Archive each analysis's inputs, retained frames, configuration, network maps,
quality reports, and software version. The pipeline does not implement a visit
database, automatically partition baseline/follow-up, or fit longitudinal models.

## Keep acquisition and processing comparable

- Document scanner, acquisition, echo times, phase encoding, fieldmap pairing,
  session dates, and intervention timing outside the public repository.
- Keep denoising, censoring, surface registration, network priors, graph densities,
  and model parameters consistent where scientifically appropriate. Record any
  intentional change; a software/configuration change can resemble biological change.
- Compare retained data duration, motion, anatomical alignment, spatial coverage,
  NSI, and network-specific quality across analysis epochs.
- Use reliability checks appropriate to your dataset; agreement across Infomap
  densities, RidgePFM weights, homogeneity diagnostics, and NSI predictions
  answer different questions and are not substitutes for repeat measurements.

## Configure the analysis boundaries explicitly

`START_SESSION` is a processing/resume control, not a complete visit-selection
language. Do not assume it defines a baseline/follow-up contrast.

For a separately prepared epoch-specific CIFTI, `PFM_INPUT_CIFTI` can specify
the mapping input and `PFM_OUTDIR` a separate output root. Supply a compatible,
properly censored CIFTI, and check the module's input preparation and provenance.
Define and record the run list and censor mask used to build that input. For NSI,
keep its analyzed input aligned with the actual mapping data rather than allowing
the default all-run concatenation to stand in for an epoch-specific dataset.

`FUNC_DIRNAME="All"` loops over task folders with shared anatomy; it is not a
longitudinal model. Standard QC summary plots follow scan order and may include
multiple runs from only one session. Never present that alone as evidence of
change over visits.

## Use individual maps for stimulation research

Run TANS on a reviewed map, with an explicit cortical search mask, coil model,
and hardware-matched calibration. Use a separate subject-relative
`TANS_OUTPUT_DIRNAME` for each distinct map/epoch, for example `tans_baseline`
and `tans_followup`. Keep each operator package intact with its target/pose
provenance; save actual session records separately from the immutable package.

Observed rMT and delivered intensity are separate observations. A planning
preview is neither an observed threshold nor a treatment record. TANS models
electric-field exposure and network overlap; it does not by itself establish
clinical suitability, safety, a therapeutic dose, or intervention efficacy.

See [TANS integration](tans-integration.md), [quality control](qc-modules.md), and
[network-output interpretation](pfm-outputs.md).

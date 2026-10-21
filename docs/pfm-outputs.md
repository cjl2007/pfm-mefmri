# PFM Outputs

[See example maps and continuous representations](output-gallery.md#individual-network-representations).
Network labels and weights can feed the [optional TANS module](tans-integration.md)
after QC and explicit selection of the input maps and metadata.

PFM runs after CIFTI concatenation and FD censoring. The default input is the
final concatenated dtseries:

```text
func/<FUNC_DIRNAME>/ConcatenatedCiftis/<FuncPrefix>_<InputTag>_Concatenated+FDlt<FD>.dtseries.nii
```

FD follows the historical `concatenate_scans.m` calculation independently for
each run. For **TR < 2 seconds**, the six motion parameters receive a zero-phase,
10th-order Butterworth bandstop from **0.2 Hz to Nyquist minus 0.019 Hz** before
FD is calculated. Runs with TR >= 2 seconds use unfiltered motion. FD uses
`round(2.5 / TR)` backward differences with MATLAB's rounding, a 50 mm rotation
radius, and zero FD for the leading lag frames of each run. Censoring retains
frames with **FD < CONCAT_FD_THRESHOLD** (default 0.3 mm).

The Python filter uses second-order sections for numerical stability and
MATLAB-compatible edge padding. This preserves the historical filter design;
it does not reproduce numerical artifacts of the old direct-form filter,
which can be unstable near TR = 2 seconds. Filtering requires more than 60
frames per run when TR < 2 seconds and fails explicitly for shorter runs.
This filtering applies to motion parameters, not the BOLD signal.

`FD.txt` stores the FD used for censoring at full precision. `ScanInfo.json`
records each run's TR, filter settings, and lag, and `PFMPrepSummary.json`
records the FD method. Motion grayplot QC uses the same calculation; keep its
threshold equal to `CONCAT_FD_THRESHOLD` to compare retained-frame counts.

Existing derivatives are not changed by updating the code. To adopt this FD
calculation for an already processed subject, regenerate concatenation and
the downstream analyses that consume its censored CIFTI, and refresh motion QC.

PFM writes strategy-specific outputs under:

```text
func/<FUNC_DIRNAME>/PFM/
  RidgePFM/
  Infomap/
```

Set `PFM_STRATEGY="all"` to run both RidgePFM and Infomap in one PFM stage.
Each strategy writes to its own subfolder under `PFM/`.

## RidgePFM

### Taxonomies, confidence, stripes, and volume maps

```bash
PFM_RF_TAXONOMY="Lynch2024"  # Lynch2024 | Du2024 | Hermosillo2024 | Kong2019
PFM_PRIORS_MAT=""            # empty: use the selected bundled taxonomy
PFM_SUBCORT_PRIORS_NII=""
PFM_RF_CONFIDENCE_ENABLE=1   # optional; default 0
PFM_RF_CONFIDENCE_SCAN_IDX="" # auto-select aligned concatenation metadata
PFM_RF_CONFIDENCE_FC_THRESHOLD=0.7
PFM_RF_VOLUME_ENABLE=1       # optional; default 0
PFM_RF_VOLUME_ANATOMICAL_TARGET="acpc" # acpc | native | custom
PFM_RF_VOLUME_REFERENCE=""  # target grid; native defaults to anat/T1w/T1w.nii.gz
PFM_RF_VOLUME_ACPC_TO_TARGET_XFM="" # native auto-selects xfms/acpc2orig.mat
```

Ribbon mapping always occurs in the subject's ACPC anatomical space, where the
white and pial surfaces are defined. The `native` target then applies the saved
ACPC-to-original transform and writes on the original `T1w.nii.gz` grid. A
`custom` target requires both `PFM_RF_VOLUME_REFERENCE` and
`PFM_RF_VOLUME_ACPC_TO_TARGET_XFM`. WTA output is a Workbench label volume with
the taxonomy names and colors embedded; score and confidence outputs remain
scalar 4D NIfTI files.

The four prior sets live in `res0urces/PriorsPFM/<taxonomy>`. Lynch2024 and
Hermosillo2024 supply cortical and subcortical spatial priors; Du2024 and Kong2019
supply cortical priors only. Targets without spatial support use FC alone:
no other taxonomy's subcortical labels are substituted. The availability map
identifies informative spatial-prior rows. Custom cortical priors override the
taxonomy selection. To add subcortical spatial support with custom priors, supply
an explicitly matching custom volume prior; otherwise volume assignments use FC.
Use a different `PFM_OUTDIR` for each taxonomy when retaining multiple analyses.

Confidence testing maps each eligible session independently, combining its runs.
Ten eligible sessions produce ten assignment maps. `ScanIdx.mat` must contain
one `[session, run]` row per input frame. New concatenation writes both the full
`ScanIdx.mat` and frame-aligned `ScanIdx_Censored.mat`; a legacy or custom CIFTI
must be supplied with matching metadata. Sessions with fewer than
`PFM_RF_REPEAT_MIN_FRAMES` (30) frames are excluded and recorded. At least four
eligible sessions are required by the combined confidence/state workflow.

The `Confidence/` directory contains:

| Output suffix | Meaning |
| --- | --- |
| `_SessionAssignmentFrequency.dscalar.nii` | Per-network fraction of valid session assignments, 0–1 |
| `_SessionAssignments.dlabel.nii` | Each session's WTA map |
| `_FCReproducibility.dscalar.nii` | FC-profile correlation for each of three pooled session splits |
| `_SplitAssignments.dlabel.nii` | The six pooled split assignments |
| `_MinimumValidSplitSessions.dscalar.nii` | Minimum number of valid sessions across split sides |
| `_AssignmentStates.dlabel.nii` | Four FC reproducibility × label stability states; zero is not estimable |
| `_SplitLabelMembership.dscalar.nii` | Every network observed across the six split maps |
| `_NetworkStripes.dlabel.nii` | Separate display map; stripes only for reproducible FC with competing labels |

The default split plan uses three seeded, disjoint-within-pair session halves.
The standalone confidence CLI also accepts an explicit `--split-plan`, including
the validation study's three held-session folds. Split pairs reuse sessions and
are not independent observations. Reproducible FC means the minimum of the three
full cortical profile correlations is at least 0.7, after the same cortical local
edge exclusion used for mapping. Same label means one nonzero label across all
six maps. Each split side must also have at least two valid sessions at that
grayordinate. This gives reproducible/same, reproducible/competing, lower/same, and
lower/competing states. Only the second state receives stripes. Numeric membership
maps retain alternatives that may not appear in a particular stripe. Stripe
widths are display geometry, not assignment probabilities.

Session frequencies measure empirical assignment consistency at the available
session duration. They are separate from `TemplateFit_R2` and from the fused-score
softmax; none is a calibrated probability of biological network membership.
The full-data WTA remains the main output and is never replaced by stripes.
The older run-level repeat and experimental areal-mixture workflows remain
separate opt-in controls.

Volume export writes WTA labels, fused-score weights, and (when enabled) session
frequencies under `Volume/`. Workbench maps cortex using the subject's midthickness,
white, and pial surfaces with `-ribbon-constrained`. Labels use categorical mapping;
scalar maps use ribbon weights normalized by mapped valid-surface coverage, so a
constant probability remains constant within covered cortex. Original CIFTI volume
components are resampled onto the requested grid and take precedence at overlap.
Coverage NIfTIs distinguish unmapped voxels (0), cortical ribbon (1), and original
volume targets (2). JSON sidecars preserve map names, label tables, and methods.
All surfaces, the reference, and CIFTI volume must already share physical space;
this exporter does not register between native and MNI spaces.

Current subcortical regression derives neighborhoods from zero-based CIFTI voxels
and subject surfaces, avoiding legacy distance matrices with a one-voxel origin
shift. The historical rounded-distance boundary is preserved as physical distance
strictly below `distance_mm + 0.5`. When aligned session/run metadata is available,
regression fits each run separately. The default censored input uses retained
frames for fitting. Updating code does not correct already regressed derivatives;
regenerate from pre-regression data or use checked corrected inputs.


RidgePFM is the current name for the method previously called Ridge-Fusion.
New runs use `PFM_STRATEGY="ridgepfm"`, the `PFM/RidgePFM/` folder, and
`RidgePFM_*` filenames. The strategy aliases `ridge_fusion`, `RidgeFusion`, and
`ridge-fusion` remain accepted; explicit `PFM_RF_*` settings still take precedence.
Existing result folders are not moved or renamed. The default FC weight is
**1.0** and the default spatial-prior weight is **0.05**.

Enable with:

```bash
PFM_ENABLE=1
PFM_STRATEGY="ridgepfm"
```

RidgePFM estimates network assignments from a regularized fusion of subject
functional connectivity evidence and spatial/network priors. The design goal is
to provide a faster, prior-guided PFM path that remains sensitive to
subject-specific deviations.

Important controls:

```bash
PFM_RF_OUTFILE="RidgePFM_VTX"
PFM_RF_FC_WEIGHT=1.0
PFM_RF_POSITIVE_FC_ONLY=0
PFM_RF_SPATIAL_WEIGHT=0.05
PFM_RF_SPATIAL_PRIOR_TRANSFORM="logit"
PFM_RF_LAMBDA=1
PFM_RF_LOCAL_EXCLUSION_MM=10
PFM_RF_FC_DEMEAN=0
PFM_RF_SUBCORT_REGRESS_ENABLE=1
PFM_RF_SMOOTHING_KERNEL=1.7
```

For inputs without GSR/MGTR-style global signal control, consider:

```bash
PFM_RF_FC_DEMEAN=1
```

`PFM_RF_SPATIAL_PRIOR_TRANSFORM="logit"` matches the legacy MATLAB
RidgePFM V4/V5 fusion rule. `PFM_RF_POSITIVE_FC_ONLY=1` is available as an
experimental V5-style diagnostic mode that fits the ridge model only on positive
target FC edges.

Primary outputs:

```text
PFM/RidgePFM/RidgePFM_VTX.dlabel.nii
PFM/RidgePFM/RidgePFM_VTX.L.border
PFM/RidgePFM/RidgePFM_VTX.R.border
PFM/RidgePFM/RidgePFM_VTX_TemplateFit_R2.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX_TopFusedScoreWeight.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX_SpatialPriorAvailable.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX_FusedScoreSoftmax.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.dlabel.nii
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.L.border
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.R.border
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation_FusedScoreSoftmax.dtseries.nii
PFM/RidgePFM/RidgePFM_Vertexwise_Review.scene
PFM/RidgePFM/RidgePFM_ArealParcellation_Review.scene
PFM/RidgePFM/RepeatAssignment/RidgePFM_VTX_PerScan_WTA.dlabel.nii
PFM/RidgePFM/RepeatAssignment/RidgePFM_VTX_ScanRepeatNetworkFrequencies.dscalar.nii
PFM/RidgePFM/RepeatAssignment/RidgePFM_VTX_ScanRepeatSupport.dscalar.nii
PFM/RidgePFM/RepeatAssignment/RidgePFM_VTX_ScanRepeatCounts.npz
PFM/RidgePFM/RepeatAssignment/run_inventory.json
PFM/RidgePFM/RepeatAssignment/status.json
PFM/RidgePFM/prep/
PFM/RidgePFM/HomogeneityTest/
```

Open `RidgePFM_VTX.dlabel.nii` in Connectome Workbench with the subject
surfaces. Template-fit R2, spatial-prior availability, and normalized fused-score
weights are separate supporting layers. Fused-score softmax values rank competing
labels but are not calibrated probabilities of true network identity or temporal
stability. `RidgePFM_Vertexwise_Review.scene` loads the vertexwise WTA/score pair
with vertexwise black border lines only. `RidgePFM_ArealParcellation_Review.scene`
loads the areal WTA/score pair with areal black border lines only.

With `PFM_RF_REPEAT_SUPPORT_ENABLE=1` and at least two runs retaining
`PFM_RF_REPEAT_MIN_FRAMES` frames, RidgePFM also
maps each run independently and averages their categorical assignments. The full
concatenated map remains the primary estimate. `ScanRepeatSupport` distinguishes
support for that primary label, competing-label frequency, modal-label support,
unavailable-run frequency, available-run count, and normalized label entropy.
`ScanRepeatNetworkFrequencies` contains the conditional frequency of every
network among available runs. These are observed temporal repeat frequencies,
not probabilities of biological membership, simultaneous membership, or a
calibrated forecast for a future scan. With fewer than two eligible runs,
`status.json` records `not_estimable` and no confidence map is invented.

Figure slots:

```text
docs/assets/pfm-ridge-fusion-dlabel.png
docs/assets/pfm-ridge-fusion-probmaps.png
```

## Infomap

Enable with:

```bash
PFM_ENABLE=1
PFM_STRATEGY="infomap"
```

Infomap preserves the community-detection workflow used by the group, including
graph-density sweeps and optional community-to-network mapping.

Important controls:

```bash
PFM_INFOMAP_GRAPH_DENSITIES_EXPR="0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001"
PFM_INFOMAP_NUM_REPS_EXPR="1,2,5,10,20,30,50,75,100"
PFM_INFOMAP_MIN_DISTANCE=30
PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1
PFM_INFOMAP_LABEL_OUTFILE="InfomapNetworkLabels"
```

If `infomap` is not on `PATH`, set:

```bash
PFM_INFOMAP_BINARY="/full/path/to/infomap"
```

For wiring tests without heavy computation:

```bash
PFM_INFOMAP_DRY_RUN=1
```

Core outputs:

```text
PFM/Infomap/Bipartite_PhysicalCommunities.dtseries.nii
PFM/Infomap/InfomapNetworkLabels_ModeConsensus.dlabel.nii
PFM/Infomap/InfomapNetworkLabels_ProbabilityConsensus.dscalar.nii
PFM/Infomap/InfomapNetworkLabels_ManualCorrections.csv
PFM/Infomap/InfomapNetworkLabels_ManualCorrections.xlsx
PFM/Infomap/HomogeneityTest/
```

Per-density outputs:

```text
PFM/Infomap/GraphDensity_<value>/
  InfomapNetworkLabels_Density<value>.dlabel.nii
  InfomapNetworkLabels_Density<value>_Confidence.dscalar.nii
  InfomapNetworkLabels_Density<value>_CommunityTable.csv
  InfomapNetworkLabels_Density<value>_AmbiguousCommunities.csv
  Bipartite_PhysicalCommunities+AlgorithmicLabeling.dlabel.nii
  Bipartite_PhysicalCommunities+AlgorithmicLabeling.L.border
  Bipartite_PhysicalCommunities+AlgorithmicLabeling.R.border
  Bipartite_PhysicalCommunities+AlgorithmicLabeling_InfoMapCommunities.dlabel.nii
  Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_WholeBrain.dtseries.nii
  Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_btwn_InfoMapCommunities.dtseries.nii
  Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv
  Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualReview.scene
```

Review priorities:

- Inspect density-specific community maps before relying on consensus.
- Use confidence and ambiguous-community tables to identify assignments that
  need manual review.
- Compare FC and spatial evidence when a community maps poorly to canonical
  network IDs.
- Use homogeneity diagnostics to help choose or justify graph density.

Per-density manual-review scene generation is controlled by:

```bash
PFM_INFOMAP_REVIEW_SCENE_ENABLE=1
PFM_INFOMAP_REVIEW_SCENE_TEMPLATE="$MEDIR/res0urces/InfomapReviewScene/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualReview.template.scene"
```

When a template is available, the labeler writes
`Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualReview.scene` into each
`GraphDensity_<value>/` folder. The scene gathers algorithmic labels, raw
Infomap community labels, whole-brain community FC, and community-to-community
FC in one Workbench review workspace.

Scene generation requires `PFM_INFOMAP_NETWORK_MAPPING_ENABLE=1`,
`PFM_INFOMAP_LABEL_WRITE_FC=1`, and a non-dry-run Infomap labeling pass because
those settings create the review files referenced by the scene.

Figure slots:

```text
docs/assets/pfm-infomap-density-map.png
docs/assets/pfm-infomap-mode-consensus.png
docs/assets/pfm-infomap-confidence.png
docs/assets/pfm-infomap-review-table.png
```

## Manual Infomap Updates

Manual edits can be applied without rerunning Infomap:

```bash
PFM_INFOMAP_UPDATE_ENABLE=1
START_FROM_MODULE="pfm_update"
```

The update module scans:

```bash
PFM_INFOMAP_UPDATE_TABLE_GLOB="GraphDensity_*/Bipartite_PhysicalCommunities+AlgorithmicLabeling_ManualCorrections.csv"
```

and writes adjusted density and consensus labels, for example:

```text
PFM/Infomap/InfomapNetworkLabels_ManualAdjusted_ModeConsensus.dlabel.nii
PFM/Infomap/InfomapNetworkLabels_ManualAdjusted_ProbabilityConsensus.dscalar.nii
```

## Areal Parcellation

Enable with:

```bash
PFM_AREAL_ENABLE=1
```

Areal parcellation sub-parcellates the network-level map into smaller areal
parcels. The output name is derived from the network-label prefix unless
overridden:

```text
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.dlabel.nii
PFM/Infomap/InfomapNetworkLabels+ArealParcellation.dlabel.nii
```

Controls:

```bash
PFM_AREAL_MIN_SIZE=30
PFM_AREAL_OUTFILE=""
```

Figure slot:

```text
docs/assets/pfm-areal-parcellation.png
```

## Experimental Multi-Network RidgePFM Maps

The opt-in multi-network extension asks whether a fixed RidgePFM parcel is
described better by two or three network FC templates than by its WTA network
alone. Enable it only with RidgePFM areal parcellation:

```bash
PFM_STRATEGY="ridgepfm"
PFM_AREAL_ENABLE=1
PFM_RF_MULTINETWORK_EXPERIMENTAL_ENABLE=1
```

Preliminary internal evaluations are encouraging, but the method still needs
broader testing across participants and acquisition conditions. It therefore
remains disabled by default and should currently be interpreted as an
experimental descriptive layer.

The experiment randomly distributes contiguous 25-frame blocks from every
scan across thirds, selects a convex one-, two-, or three-network description
on two thirds, and scores it on the held-out third. The default eight repeats
produce 24 held-out folds. Added networks must recur in at least 80% of folds,
improve held-out fit in at least 80%, have positive mean gain in every split
repeat, and carry a nontrivial mixture weight (defaults: 0.20 for the second
network and 0.15 for the third).

The default parcelwise null preserves each parcel's complete held-out gain
landscape but permutes the correspondence between the added-network identity
selected on training data and its held-out identity. Thus generic improvement
from simply adding model flexibility is included in the null. Two thousand
repeat-preserving permutations yield parcelwise p-values, followed by
Benjamini-Hochberg FDR at 0.05 separately in cortex and subcortex.
Three-network support must pass both sequential comparisons (`2 vs 1` and `3
vs 2`). The older whole-brain max-statistic result is retained in the CSV and
dscalar as a separate extreme-confidence FWER tier.

This layer keeps the full-data parcel geometry fixed and does not rerun areal
parcellation or RidgePFM within each fold. It is a held-out FC-template
mixture diagnostic, so it supplements rather than replaces the original WTA
and probability maps.

Important controls:

```bash
PFM_RF_MULTINETWORK_STRIPE_PRESET="balanced" # strict | balanced | loose | custom
```

This is the recommended single tightness knob:

| Preset | Parcel decision | Mixture/stability | Local stripe extent |
| --- | --- | --- | --- |
| `strict` | Whole-brain FWER tier | 90% folds; weights 0.25/0.20 | `p>=0.20`, ratio `>=0.667`, clusters 20/5 |
| `balanced` | Domain-specific FDR 0.05 | 80% folds; weights 0.20/0.15 | `p>=0.10`, ratio `>=0.40`, clusters 20/5 |
| `loose` | Domain-specific FDR 0.10 | 70% folds; weights 0.15/0.10 | `p>=0.075`, ratio `>=0.333`, clusters 10/3 |
| `custom` | Uses the individual settings below | User-defined | User-defined |

`balanced` is the recommended experimental starting point. The preset is an
explicit scientific operating-point choice, not a cosmetic stripe-width
setting. The null/FDR calculations remain data-driven within that choice, and
the best operating point should be reassessed as broader validation accrues.

Advanced controls. Split/null-generation settings remain active for every
preset; `custom` additionally uses the individual tightness thresholds:

```bash
PFM_RF_MULTINETWORK_SPLIT_REPEATS=8
PFM_RF_MULTINETWORK_SPLIT_BLOCK_LENGTH=25
PFM_RF_MULTINETWORK_MIN_CONSISTENCY_FRACTION=0.8
PFM_RF_MULTINETWORK_MIN_POSITIVE_GAIN_FRACTION=0.8
PFM_RF_MULTINETWORK_SELECTION_CORRECTION="fdr"
PFM_RF_MULTINETWORK_FDR_ALPHA=0.05
PFM_RF_MULTINETWORK_PARCEL_NULL_ITERATIONS=2000
PFM_RF_MULTINETWORK_MIN_SECOND_WEIGHT=0.20
PFM_RF_MULTINETWORK_MIN_THIRD_WEIGHT=0.15
PFM_RF_MULTINETWORK_NULL_ALPHA=0.05
PFM_RF_MULTINETWORK_MAX_NETWORKS=3
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_PROB_MIN=0.10
PFM_RF_MULTINETWORK_LOCAL_SECONDARY_RATIO_MIN=0.40
PFM_RF_MULTINETWORK_MIN_CORTICAL_CLUSTER_VERTICES=20
PFM_RF_MULTINETWORK_MIN_SUBCORTICAL_CLUSTER_VOXELS=5
```

Outputs are written under `PFM/RidgePFM/ExperimentalMultiNetwork/`:

```text
RidgePFM_VTX+ExperimentalMultiNetworkCV_ParcelSummary.csv
RidgePFM_VTX+ExperimentalMultiNetworkCV_FoldResults.csv
RidgePFM_VTX+ExperimentalMultiNetworkCV_Diagnostics.dscalar.nii
RidgePFM_VTX+ExperimentalMultiNetworkCV_Summary.json
RidgePFM_VTX+ExperimentalMultiNetwork_Parcel_NetworkStripes.dlabel.nii
RidgePFM_VTX+ExperimentalMultiNetwork_Parcel_SubcorticalNetworkStripes_0p5mm.nii.gz
RidgePFM_VTX+ExperimentalMultiNetwork_Localized_NetworkStripes.dlabel.nii
RidgePFM_VTX+ExperimentalMultiNetwork_Localized_Diagnostics.dscalar.nii
RidgePFM_VTX+ExperimentalMultiNetwork_Localized_ClusterSummary.csv
```

The parcel dlabel stripes every null-selected parcel that also meets the
configured cortical/subcortical parcel-size requirement. The localized dlabel
is more conservative spatially: within a selected parcel, it requires the same
network set to have local RidgePFM support and form a minimum-size
contiguous surface or volume cluster. Thus a parcel may be striped in the
parcel map but remain WTA in the localized map when its evidence is diffuse.
This is a parcel-anchored vertex/voxel translation, not a separate vertexwise
CV/null refit. Outside retained local clusters it exactly preserves the original
RidgePFM vertex/voxel WTA labels.
The 0.10 local cutoff is twice the uniform 1/20 reference probability; it is an
effect-size screen, not a calibrated posterior-probability threshold.

The candidate-identity null is conditional on the observed candidate-gain
landscape; it is not a generative residual surrogate under a fitted
one-network timeseries model. Treat this as an experimental descriptive layer
until it replicates across subjects and a full residual-surrogate or per-fold
RidgePFM analysis confirms calibration.

## Workbench Review Checklist

For each subject:

- Open the subject midthickness/inflated surfaces.
- Load the network dlabel, borders, R2/probability maps, and homogeneity figure.
- Check whether maps respect expected sensory/motor/association topography.
- Check for session-specific artifacts using QC reports before interpreting
  unusual network deviations.
- Compare RidgePFM and Infomap when both were run.
- Record any manual Infomap corrections in the per-density correction tables.

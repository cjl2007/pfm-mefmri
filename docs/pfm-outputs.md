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
PFM_RF_LAMBDA=10
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
PFM/RidgePFM/RidgePFM_VTX_R2.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX_ProbMaps.dtseries.nii
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.dlabel.nii
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.L.border
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation.R.border
PFM/RidgePFM/RidgePFM_VTX+ArealParcellation_ProbMaps.dtseries.nii
PFM/RidgePFM/RidgePFM_Vertexwise_Review.scene
PFM/RidgePFM/RidgePFM_ArealParcellation_Review.scene
PFM/RidgePFM/prep/
PFM/RidgePFM/HomogeneityTest/
```

Open `RidgePFM_VTX.dlabel.nii` in Connectome Workbench with the subject
surfaces. Use `RidgePFM_VTX_R2.dtseries.nii` and probability maps as
confidence-style supporting layers. `RidgePFM_Vertexwise_Review.scene` loads
the vertexwise WTA/probability pair with vertexwise black border lines only.
`RidgePFM_ArealParcellation_Review.scene` loads the areal WTA/probability
pair with areal black border lines only.

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

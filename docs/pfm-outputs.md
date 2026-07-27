# PFM Outputs

PFM runs after CIFTI concatenation and FD censoring. The default input is the
final concatenated dtseries:

```text
func/<FUNC_DIRNAME>/ConcatenatedCiftis/<FuncPrefix>_<InputTag>_Concatenated+FDlt<FD>.dtseries.nii
```

PFM writes strategy-specific outputs under:

```text
func/<FUNC_DIRNAME>/PFM/
  RidgeFusion/
  Infomap/
```

Set `PFM_STRATEGY="all"` to run both Ridge-Fusion and Infomap in one PFM stage.
Each strategy writes to its own subfolder under `PFM/`.

## Ridge-Fusion

Enable with:

```bash
PFM_ENABLE=1
PFM_STRATEGY="ridge_fusion"
```

Ridge-Fusion estimates network assignments from a regularized fusion of subject
functional connectivity evidence and spatial/network priors. The design goal is
to provide a faster, prior-guided PFM path that remains sensitive to
subject-specific deviations.

Important controls:

```bash
PFM_RF_OUTFILE="RidgeFusion_VTX"
PFM_RF_FC_WEIGHT=1.0
PFM_RF_SPATIAL_WEIGHT=0.1
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

Primary outputs:

```text
PFM/RidgeFusion/RidgeFusion_VTX.dlabel.nii
PFM/RidgeFusion/RidgeFusion_VTX.L.border
PFM/RidgeFusion/RidgeFusion_VTX.R.border
PFM/RidgeFusion/RidgeFusion_VTX_R2.dtseries.nii
PFM/RidgeFusion/RidgeFusion_VTX_ProbMaps.dtseries.nii
PFM/RidgeFusion/prep/
PFM/RidgeFusion/HomogeneityTest/
```

Open `RidgeFusion_VTX.dlabel.nii` in Connectome Workbench with the subject
surfaces. Use `RidgeFusion_VTX_R2.dtseries.nii` and probability maps as
confidence-style supporting layers.

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
  InfomapNetworkLabels_Density<value>_FC.dscalar.nii
  InfomapNetworkLabels_Density<value>_CommunityTable.csv
  InfomapNetworkLabels_Density<value>_AmbiguousCommunities.csv
  Bipartite_PhysicalCommunities+AlgorithmicLabeling.dlabel.nii
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
PFM/RidgeFusion/RidgeFusion_VTX+ArealParcellation.dlabel.nii
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

## Experimental Multi-Network Ridge-Fusion Maps

The opt-in multi-network extension asks whether a fixed Ridge-Fusion parcel is
described better by two or three network FC templates than by its WTA network
alone. Enable it only with Ridge-Fusion areal parcellation:

```bash
PFM_STRATEGY="ridge_fusion"
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
parcellation or Ridge-Fusion within each fold. It is a held-out FC-template
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

Outputs are written under `PFM/RidgeFusion/ExperimentalMultiNetwork/`:

```text
RidgeFusion_VTX+ExperimentalMultiNetworkCV_ParcelSummary.csv
RidgeFusion_VTX+ExperimentalMultiNetworkCV_FoldResults.csv
RidgeFusion_VTX+ExperimentalMultiNetworkCV_Diagnostics.dscalar.nii
RidgeFusion_VTX+ExperimentalMultiNetworkCV_Summary.json
RidgeFusion_VTX+ExperimentalMultiNetwork_Parcel_NetworkStripes.dlabel.nii
RidgeFusion_VTX+ExperimentalMultiNetwork_Parcel_SubcorticalNetworkStripes_0p5mm.nii.gz
RidgeFusion_VTX+ExperimentalMultiNetwork_Localized_NetworkStripes.dlabel.nii
RidgeFusion_VTX+ExperimentalMultiNetwork_Localized_Diagnostics.dscalar.nii
RidgeFusion_VTX+ExperimentalMultiNetwork_Localized_ClusterSummary.csv
```

The parcel dlabel stripes every null-selected parcel that also meets the
configured cortical/subcortical parcel-size requirement. The localized dlabel
is more conservative spatially: within a selected parcel, it requires the same
network set to have local Ridge-Fusion support and form a minimum-size
contiguous surface or volume cluster. Thus a parcel may be striped in the
parcel map but remain WTA in the localized map when its evidence is diffuse.
This is a parcel-anchored vertex/voxel translation, not a separate vertexwise
CV/null refit. Outside retained local clusters it exactly preserves the original
Ridge-Fusion vertex/voxel WTA labels.
The 0.10 local cutoff is twice the uniform 1/20 reference probability; it is an
effect-size screen, not a calibrated posterior-probability threshold.

The candidate-identity null is conditional on the observed candidate-gain
landscape; it is not a generative residual surrogate under a fitted
one-network timeseries model. Treat this as an experimental descriptive layer
until it replicates across subjects and a full residual-surrogate or per-fold
Ridge-Fusion analysis confirms calibration.

## Workbench Review Checklist

For each subject:

- Open the subject midthickness/inflated surfaces.
- Load the network dlabel, borders, R2/probability maps, and homogeneity figure.
- Check whether maps respect expected sensory/motor/association topography.
- Check for session-specific artifacts using QC reports before interpreting
  unusual network deviations.
- Compare Ridge-Fusion and Infomap when both were run.
- Record any manual Infomap corrections in the per-density correction tables.

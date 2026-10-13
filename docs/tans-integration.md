# TANS as the final pfm-mefmri module

PFM-mefmri owns acquisition processing, QC, and network mapping. The separate
[pfm-tans repository](https://github.com/cjl2007/pfm-tans) owns head-model
preparation, target preparation, coil search, dose reports, and operator exports.
`modules/mefmri_func_tans.sh` is the connection, not a second TANS implementation.

## Install and configure

1. Install pfm-tans and its documented SimNIBS/Workbench/Python dependencies.
   Keep its `pipeline/`, `search/`, and `res0urces/` directories together.
2. First create a complete [pfm-mefmri study config](configuration.md).
   Copy [tans_study.example.yaml](../config/tans_study.example.yaml) into a
   study-owned location. Set software paths, input maps, coils, calibration
   identities, search settings, retention, robustness, and dose limits.
3. Add [tans_module.example.sh](../config/tans_module.example.sh) to the end of
   a complete pfm-mefmri study config, by copying its settings or sourcing a
   study-owned copy. It is a configuration fragment, not a launcher.
4. Set `TANS_CODE_ROOT` to `/path/to/pfm-tans/pipeline`, **including `pipeline`**,
   and `TANS_CONFIG_FILE` to the edited YAML. `TANS_PYTHON` must provide the
   dependencies required by the TANS driver; the SimNIBS environment is the
   documented example.

The release and fragment templates keep `TANS_ENABLE=0`. Replace every
`/path/to` and `REPLACE_WITH` placeholder in the YAML, including both matched
hardware calibration IDs. No coil/stimulator pair is selected for you. The
YAML `simulation.coil` field should match `simulation.treatment_coil`.
The default handle mode is `as_simulated`; only enable automatic 180-degree
flipping after explicitly verifying the coil's magnitude equivalence.

For example, after configuring the YAML, set:

```bash
PFM_ENABLE=1  # for a full upstream run; requires completed PFM for a TANS-only restart
TANS_ENABLE=1
TANS_CODE_ROOT="/opt/pfm-tans/pipeline"
TANS_CONFIG_FILE="/data/MyStudy/config/tans.yaml"
TANS_PYTHON="/opt/SimNIBS-4.5/simnibs_env/bin/python"
TANS_TARGET_NETWORK="CinguloOpercular"  # illustrative target in the bundled metadata
TANS_OUTPUT_DIRNAME="tans"
STOP_AFTER_MODULE="tans"
```

For a full processing chain, keep the appropriate upstream start stage and
enable PFM. For completed anatomy and PFM, also use:

```bash
START_FROM_MODULE="tans"
```

Launch through the normal pipeline entry point:

```bash
bash /opt/pfm-mefmri/bin/mefmri_pipeline.sh \
  /data/MyStudy/SUB001 /data/MyStudy/config/mefmri.sh
```

The shell `TANS_*` settings passed by the wrapper override their corresponding
YAML settings. Empty coil-path and calibration-ID overrides leave hardware
selection to the study YAML and TANS's exact-match resolution. Explicitly supplied
values still take precedence. Never rely on example coils or calibration tables
for different physical hardware, pulse shapes, or current directions.

## Inputs and map selection

The wrapper checks native `anat/T1w/T1w.nii.gz`, native `anat/T2w/T2w.nii.gz`,
the YAML, driver, and PFM folder. TANS additionally requires compatible HCP-style
fs_LR 32k surfaces/areas, anatomical labels and native/ACPC transforms. See its
full input contract. Passing the wrapper's checks alone is not input validation.

Default mapping input root:

```text
<SubjectDir>/func/<FUNC_DIRNAME>/PFM/RidgePFM/
```

The example YAML expects `RidgePFM_VTX.dlabel.nii` and
`RidgePFM_VTX_ProbMaps.dtseries.nii`. For Infomap:

- Set `TANS_PFM_DIR` to the chosen participant's `PFM/Infomap` directory.
- Set YAML `inputs.categorical_networks` to
  `"{pfm_dir}/InfomapNetworkLabels_ModeConsensus.dlabel.nii"`.
- Set YAML `inputs.network_probabilities` to
  `"{pfm_dir}/InfomapNetworkLabels_ProbabilityConsensus.dscalar.nii"`.
- Use the manual-adjusted counterparts when that is the reviewed result.
  Confirm that map order and network metadata agree; unlabeled raw Infomap
  community IDs are not automatically canonical network identities.

For existing `PFM/RidgeFusion/` results, the wrapper falls back to that folder
when `PFM/RidgePFM/` is absent. Keep the legacy `RidgeFusion_*` filenames in
your study YAML, or supply explicit paths to the reviewed maps. If both folders
exist, the default is `RidgePFM`; use `TANS_PFM_DIR` to select another source.

Changing `PFM_STRATEGY` does not automatically select TANS filenames or metadata.
Custom WTA/weighted activation or FC maps need an explicit YAML input contract.
WTA-only inputs can supply 0/1 stand-in weights; weighted-only inputs require an
explicit percentile cutoff to derive support. These fallbacks do not recover
missing uncertainty or calibrate probabilities. A motor-reference workflow also
needs an SM-Hand target. With only a custom nonmotor target, use
`TANS_INCLUDE_SM_HAND_CALIBRATION=0`; do not imply that motor-calibrated dose
pages can be generated without a compatible motor reference.

## Study choices that deserve explicit review

| Setting | Meaning |
| --- | --- |
| `TANS_TARGET_NETWORK` | Network name/index in the supplied metadata, not an anatomical region. |
| `TANS_TARGET_SEARCH_SPACE` | Default is a left-prefrontal recipe. Supply an appropriate binary CIFTI mask for another cortical zone. |
| `TANS_INCLUDE_SM_HAND_CALIBRATION` | Include the separate motor-reference target when the required inputs exist. |
| `TANS_SM_HAND_COIL_MODEL`, `TANS_TREATMENT_COIL_MODEL` | Optional shell overrides for the actual calibration and treatment coils. |
| `TANS_SM_HAND_DIDT_CALIBRATION`, `TANS_TREATMENT_DIDT_CALIBRATION` | Exact hardware/pulse calibration identities or tables. |
| `TANS_ROBUSTNESS_ENABLE` | Evaluate placement-error robustness (`1`) or use the nominal fine-search winner (`0`). |
| `TANS_MAXIMUM_TREATMENT_RMT_RATIO` | Configured upward console-output limit; explicitly empty disables it. Not a validated safety rule. |
| `TANS_START_FROM`, `TANS_STOP_AFTER`, `TANS_RESUME` | TANS-internal stages and reuse; distinct from pfm-mefmri module routing. |
| `TANS_OUTPUT_DIRNAME` | One directory name under the participant root; use separate names for different map/epoch analyses. |

The selected network and cortical zone can be changed independently. The SAL
example is not a restriction on the algorithm. Feasibility still depends on
anatomy, access, field reach, data quality, and study protocol.

## Outputs and observed-rMT workflow

Default outputs are under `<SubjectDir>/tans/`:

```text
reports/index.html                       research/QC entry point
targets/<target>/reports/                target and search reports
shared/export_<subject>/index.html       portable operator package
shared/export_<subject>/<target>/dose-calibration/DoseCalibrationReport.html
```

TANS builds its native-input head model separately; it does not reuse the
upstream anatomy mask-generation CHARM model as if the coordinate frames were
interchangeable. Review all registration and navigation exports.

Before observed rMT is confirmed, dosing is a preview. Confirm the measured rMT
and matched motor-reference setup, inspect modeled spread, and separately record
the actual treatment intensity. Browser-generated session downloads must be
saved explicitly outside the immutable export package.

For `FUNC_DIRNAME="All"`, do not enable a shared subject-level TANS output for
multiple different task maps. Run TANS separately on the selected task/map source
with distinct output names. TANS will reject incompatible reuse rather than
silently merge different targets or map provenance.

## Verification scope

`tests/test_tans_module.py` checks quoting, paths, options, and YAML hardware
fallback. `tests/test_tans_pipeline_hook.py` checks module routing with stub
scientific stages and a stub TANS driver. Neither is a new head-model/FEM run,
an end-to-end scanner-data validation, or clinical validation.

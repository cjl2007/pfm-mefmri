# Task modeling

Experimental — work in progress.

Enable `TASK_MODEL_ENABLE=1` and set `STOP_AFTER_MODULE="pfm"` to run
per-run task GLMs after selected-dataset PFM and optional TANS.
The runner adds only the optional stage dispatch; existing preprocessing,
denoising, mapping, and PFM modules are unchanged. This module is
opt-in; existing configurations keep their previous behavior. Resume with
`PIPELINE_RUN_CORE=0` to reuse completed preprocessing and PFM.

Reusable recipes live in `config/task_models/*.yaml`. Schema 1 currently
supports the `block_events` adapter: a two-condition block design, implicit
rest, canonical SPM HRF plus time derivative, AR1 errors, motion24, raw Power
FD spike regressors, eroded WM/CSF means, and cosine drift. The motor recipe
compares right thumb flexion with index abduction. Other task designs need a
supported adapter and recipe; the effort prototype is not registered.

`TASK_MODEL_TASKS=("motor")` explicitly selects a recipe. When blank, supported
recipes are selected from the wrapper's functional tasks; rest is skipped
unless a matching recipe exists. Missing inputs for a selected supported task
are errors. Modeling never uses the PFM concatenation or MGTR data: OCME and
OCME+MEICA are fitted independently using the original per-run frame clock.

Place the run's events and timing metadata under:

```text
func/unprocessed/motor/session_1/run_1/
  Motor_S1_R1_events.csv
  Motor_S1_R1_taskmodel.json
```

CSV columns: `onset`, `duration`, `trial_type`, `phase`, `hand`. Conditions are
`thumb`, `index`, and `rest`; modeled movement blocks use `phase=experimental`.
Use observed display onsets/durations, in seconds. Practice is not included.
Task timing is relative to the first accepted trigger, with an explicit run
metadata offset to the first saved raw volume. It must not be inferred from
filenames or the number of logged triggers. Timing metadata example:

```json
{
  "alignment_status": "confirmed",
  "onset_offset_seconds": 0,
  "slice_reference": "pipeline_slice_timing_provenance",
  "source": "Acquisition log and verified first-trigger correspondence"
}
```

`assumed_authorized` is also supported, with the user's assumption and relevant
timing QC recorded. Pending alignment refuses fitting. A positive offset means
the task's trigger zero occurred after the first raw saved volume. The backend
subtracts `rmVols.txt * TR` for actual preprocessing removals, requires raw and
processed frame counts to agree, and derives the effective sampling reference from raw JSON slice times plus
the recorded fractional-TR `SliceTiming.txt` shifts. It verifies the pipeline's
rank-group convention and records the residual spread across slice groups. Other preprocessing timing conventions need an
explicit backend update. Current tissue/scene support requires T1w atlas space.

Outputs are versioned under:

```text
task_models/
  index.html
  motor/session_1/run_1/<input-and-code-fingerprint>/
    index.html
    recipe.yaml
    timing.json
    input_provenance.json
    OCME/
    OCME+MEICA/
    complete.json
```

Each input has effect, standard error, t, signed z, two-sided p, and two-sided
BH FDR q maps; design/events tables; FD and design diagnostics; residual maps;
and bilateral inflated/flatmap scenes and PNGs. The offline HTML report toggles
contrasts and unthresholded/FDR displays with matched color ranges and visible
white individual network borders. Network labels default to the selected
PFM root's `RidgePFM/RidgePFM_VTX.dlabel.nii`; set
`TASK_MODEL_NETWORK_LABEL` explicitly for another labeling. Adjacent `.L.border`
and `.R.border` files must exist. Combined rest+task PFM borders are descriptive;
they are not an independent network-specificity test.

FDR is per contrast/input over valid grayordinates, not across all contrasts.
The reverse digit contrast is the sign reversal of the primary contrast.
Differences between the two inputs' thresholded maps do not test a formal
between-input difference. Digit effects also include differences between the
two movement types. This is a single-run localization analysis, not a claim
of independent replication.

Complete outputs are reused only when their input/config/code fingerprints
and output hashes match. Existing incomplete or edited results are preserved
and cause an error; use a separate `TASK_MODEL_OUTDIR` to rerun them. No fitting
result is overwritten. The run-level directories also permit additional runs
without pretending that a pooled model has been fitted.

Backend references: [Nilearn design matrix](https://nilearn.github.io/stable/modules/generated/nilearn.glm.first_level.make_first_level_design_matrix.html)
and [Nilearn AR GLM](https://nilearn.github.io/stable/modules/generated/nilearn.glm.first_level.run_glm.html).

See [independent postprocessing controls](postprocessing.md). TANS and task modeling are enabled independently of core start/stop routing.

Task selection is a Bash string array, for example `TASK_MODEL_TASKS=("motor" "effort")`. Each entry resolves exactly to `<TASK_MODEL_CONFIG_DIR>/<task>.yaml`, whose `task` field must match. Every requested recipe must exist and have a supported adapter. Duplicate or unsafe names are errors. Arrays are serialized as JSON across process boundaries so no names are silently dropped. Only motor is currently provided as a supported recipe; the effort prototype is not executable.

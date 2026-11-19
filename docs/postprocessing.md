# Independent optional postprocessing

`START_FROM_MODULE` and `STOP_AFTER_MODULE` select the ordered core pipeline,
from `validate` through `pfm`. They do not select TANS or task modeling.

```bash
PIPELINE_RUN_CORE=1
START_FROM_MODULE="validate"
STOP_AFTER_MODULE="pfm"
TANS_ENABLE=0
TASK_MODEL_ENABLE=1
TASK_MODEL_TASKS=("motor")
TASK_MODEL_CONFIG_DIR="$MEDIR/config/task_models"
```

After successful completion of the requested core range, each enabled consumer
runs once. Enabling task modeling does not enable PFM or TANS. Each consumer
checks its own existing inputs. TANS requires the selected PFM output; task
GLMs require processed per-run inputs and, if borders are enabled, the specified
network labels and surfaces. Missing inputs are errors, not requests to rerun
preprocessing or create PFM outputs.

Run only consumers against existing data:

```bash
PIPELINE_RUN_CORE=0
TANS_ENABLE=0
TASK_MODEL_ENABLE=1
TASK_MODEL_TASKS=("motor")
TASK_MODEL_CONFIG_DIR="$MEDIR/config/task_models"
# Start/stop settings are ignored when PIPELINE_RUN_CORE=0.
```

Task names must be explicit when modeling is enabled. Each must have a matching
supported YAML in the configured directory. No automatic selection of every
supported task occurs. Per-run events and timing metadata remain in unprocessed.
The task adapter and scientific model are described in [task modeling](task-modeling.md).

The optional consumers are independent: a failing TANS invocation is recorded,
and task modeling is still attempted if enabled. Any failure makes the overall
command unsuccessful. Core failure prevents optional consumers from running.
Logs live in `<SubjectDir>/QC/Postprocessing/ModuleLogs`.

For a validation-only check, turn both consumers off; otherwise enabled consumers
will still be attempted after validation, independently of the core stop point.
Old `START_FROM_MODULE=tans` or `task_model` settings should be replaced with
`PIPELINE_RUN_CORE=0` and the appropriate enable switches. Old optional stop
values should be replaced by the desired core stop (usually `pfm`).

Task selection is a Bash string array, for example `TASK_MODEL_TASKS=("motor" "effort")`. Each entry resolves exactly to `<TASK_MODEL_CONFIG_DIR>/<task>.yaml`, whose `task` field must match. Every requested recipe must exist and have a supported adapter. Duplicate or unsafe names are errors. Arrays are serialized as JSON across process boundaries so no names are silently dropped. Only motor is currently provided as a supported recipe; the effort prototype is not executable.

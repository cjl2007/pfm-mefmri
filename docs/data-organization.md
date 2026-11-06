# Data organization

Anatomy is shared by each participant. Functional inputs and outputs are grouped
by dataset, session, and run. This example uses **rest**, **Motor**, and
**WorkingMemory**; each dataset can have a different number of acquisitions.

## Before processing

Stage your inputs with [the importer](setup-and-inputs.md#importing-inputs).
The example below uses three echoes and TOPUP fieldmaps.

```text
MyStudy/
└── SUB001/
    ├── anat/
    │   └── unprocessed/
    │       ├── T1w/
    │       │   └── T1w_1.nii.gz
    │       └── T2w/                         # optional
    │           └── T2w_1.nii.gz
    └── func/
        └── unprocessed/
            ├── rest/
            │   ├── field_maps/             # TOPUP example: matched AP/PA + JSON
            │   │   ├── AP_S1_R1.nii.gz
            │   │   ├── AP_S1_R1.json
            │   │   ├── PA_S1_R1.nii.gz
            │   │   ├── PA_S1_R1.json
            │   │   └── ...                 # additional acquired fieldmap pairs
            │   ├── session_1/
            │   │   ├── run_1/
            │   │   │   ├── Rest_S1_R1_E1.nii.gz
            │   │   │   ├── Rest_S1_R1_E1.json
            │   │   │   ├── Rest_S1_R1_E2.nii.gz
            │   │   │   ├── Rest_S1_R1_E2.json
            │   │   │   ├── Rest_S1_R1_E3.nii.gz
            │   │   │   └── Rest_S1_R1_E3.json
            │   │   └── run_2/              # Rest_S1_R2_E1...E3, each with JSON
            │   └── session_2/
            │       ├── run_1/              # Rest_S2_R1_E1...E3, each with JSON
            │       └── run_2/              # Rest_S2_R2_E1...E3, each with JSON
            ├── Motor/
            │   ├── field_maps/             # fieldmaps acquired for Motor
            │   └── session_1/
            │       ├── run_1/              # Motor_S1_R1_E1...E3, each with JSON
            │       └── run_2/              # Motor_S1_R2_E1...E3, each with JSON
            └── WorkingMemory/
                ├── field_maps/             # fieldmaps acquired for WorkingMemory
                ├── session_1/
                │   └── run_1/              # WorkingMemory_S1_R1_E1...E3 + JSON
                └── session_2/
                    └── run_1/              # WorkingMemory_S2_R1_E1...E3 + JSON
```

`E1...E3` in comments abbreviates separate E1, E2, and E3 files, each with its
own JSON sidecar. For single-echo data, supply E1 only. Folder and prefix names
are case-sensitive.

### What validate checks

- At least one `anat/unprocessed/T1w/T1w*.nii.gz`; T2w is optional.
- Numeric `session_<N>/run_<N>` folders and filenames matching the configured
  prefix, session, run, and echo indices.
- A JSON sidecar for every echo, with matching, continuous echo indices from E1.
- `EchoTime`, `RepetitionTime`, and `SliceTiming` metadata; readout and
  phase-encoding metadata when using fieldmap correction.
- Matching timepoint counts across echoes when FSL is available, plus the
  timepoint requirement for a numeric tedana PCA component count.
- Paired AP/PA fieldmaps and JSON files for the TOPUP layout shown above.

SBRefs can also be supplied; coregistration can create a fallback when absent.
Other distortion-correction strategies have different input requirements; see
[the setup guide](setup-and-inputs.md#expected-input-layout).

For these dataset names, the corresponding selection is:

```bash
FUNC_DIRNAME=("rest" "Motor" "WorkingMemory")
FUNC_FILE_PREFIX=("Rest" "Motor" "WorkingMemory")
PFM_INPUT_TASKS=("rest")
```

This processes all three datasets and selects resting-state data for PFM.
Start a fresh run from `validate` to check the staged inputs.

## After processing

Selected example outputs for an ME-ICA branch and a 0.3 mm FD threshold:

```text
SUB001/
├── anat/
│   ├── unprocessed/                        # staged anatomical inputs
│   ├── T1w/
│   │   ├── T1w_acpc_dc_restore_brain.nii.gz
│   │   └── fsaverage_LR32k/                 # individual surfaces
│   └── MNINonLinear/                       # atlas-space anatomy / transforms
├── func/
│   ├── unprocessed/                        # rest/, Motor/, WorkingMemory/
│   ├── xfms/
│   │   ├── rest/
│   │   ├── Motor/
│   │   └── WorkingMemory/
│   ├── rest/
│   │   ├── session_1/
│   │   │   ├── run_1/
│   │   │   │   ├── TE.txt
│   │   │   │   ├── SBref.nii.gz
│   │   │   │   ├── Rest_E1_acpc.nii.gz
│   │   │   │   ├── Rest_OCME+MEICA.nii.gz
│   │   │   │   └── Rest_OCME+MEICA.dtseries.nii
│   │   │   └── run_2/
│   │   ├── session_2/
│   │   │   ├── run_1/
│   │   │   └── run_2/
│   │   ├── field_maps/
│   │   ├── ConcatenatedCiftis/
│   │   │   ├── Rest_OCME+MEICA_Concatenated.dtseries.nii
│   │   │   ├── Rest_OCME+MEICA_Concatenated+FDlt0p3.dtseries.nii
│   │   │   ├── FD.txt
│   │   │   ├── ScanInfo.json
│   │   │   ├── ScanIdx.mat
│   │   │   └── ScanIdx_Censored.mat
│   │   └── QC/                             # Coreg/, Fieldmaps/, MotionGrayPlots/,
│   │                                       # NSI/, CrawlingSeedFC/, RunMetadata/
│   ├── Motor/
│   │   ├── session_1/
│   │   │   ├── run_1/                      # Motor_OCME+MEICA.dtseries.nii
│   │   │   └── run_2/
│   │   ├── field_maps/
│   │   ├── ConcatenatedCiftis/              # Motor_... + scan/frame metadata
│   │   └── QC/
│   └── WorkingMemory/
│       ├── session_1/
│       │   └── run_1/                      # WorkingMemory_OCME+MEICA.dtseries.nii
│       ├── session_2/
│       │   └── run_1/
│       ├── field_maps/
│       ├── ConcatenatedCiftis/              # WorkingMemory_... + metadata
│       └── QC/
└── pfm/                                    # selected PFM input datasets
    ├── dataset.json
    ├── ConcatenatedCiftis/
    │   ├── PFM_Combined_Concatenated+FDlt0p3.dtseries.nii
    │   ├── PFMPrepSummary.json
    │   ├── ScanInfo.json
    │   └── ScanIdx_Censored.mat
    └── RidgePFM/                           # when RidgePFM is enabled
        ├── RidgePFM_VTX.dlabel.nii
        ├── RidgePFM_VTX_FusedScoreSoftmax.dtseries.nii
        └── RidgePFM_VTX.L.border
```

Preprocessing suffixes and optional products depend on the configuration. Run
and session indices remain in the folder structure; processed filenames use
the dataset prefix. `ScanInfo.json` and `ScanIdx_Censored.mat` retain the
run/session identities associated with the concatenated and retained frames.

Each dataset has its own functional derivatives and QC. The shared subject-level
`pfm/` folder records the selected PFM inputs in `dataset.json` and stores the
enabled RidgePFM and/or Infomap results. Existing older results are not relocated.

[Back to README](../README.md#data-organization)

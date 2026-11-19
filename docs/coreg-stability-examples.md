# Within-subject alignment examples

Six individuals illustrate how refinement affects consistency across acquisitions. Each movie compares anatomy alignment with anatomy plus within-subject refinement, using the same scans and display settings.

## Differences across scans

![Scan-reference differences across six individuals](assets/coreg-stability/six_subject_stability.png)

The upper row for each group shows differences between consecutive scans. The lower row shows each scan's median difference to the others. Vertical lines mark changes of session. Lower RMS means more similar images.

<details>
<summary>Measurements</summary>

The examples use saved two-iteration, session-balanced refinement. Baseline references were checked against the saved QA images. Comparisons use identical anatomical masks, per-scan median-intensity normalization, and matched display limits. RMS includes contrast and noise effects as well as alignment.

Runs missing the required saved transforms were excluded: 14 for Individual 4 and 16 for Individual 5.

</details>

## Individual 1

![Alignment comparison for Individual 1](assets/coreg-stability/individual_1_matched_alignment.gif)

## Individual 2

![Alignment comparison for Individual 2](assets/coreg-stability/individual_2_matched_alignment.gif)

## Individual 3

![Alignment comparison for Individual 3](assets/coreg-stability/individual_3_matched_alignment.gif)

## Individual 4

![Alignment comparison for Individual 4](assets/coreg-stability/individual_4_matched_alignment.gif)

## Individual 5

![Alignment comparison for Individual 5](assets/coreg-stability/individual_5_matched_alignment.gif)

## Individual 6

![Alignment comparison for Individual 6](assets/coreg-stability/individual_6_matched_alignment.gif)

[Numerical summaries](assets/coreg-stability/subject_summaries.csv) · [Figure PDF](assets/coreg-stability/six_subject_stability.pdf) · [Longitudinal alignment](longitudinal-alignment.md)

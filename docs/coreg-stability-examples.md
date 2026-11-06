# Within-subject alignment: six-subject examples

These examples reproduce the saved two-iteration, session-balanced template-refinement demonstration, starting from native scan references and the original correction transforms. The regenerated baseline was checked against saved QA images before refinement. This does not validate every current production refinement option.

![Neighboring and median pairwise reference RMS across six participants](assets/coreg-stability/six_subject_stability.png)

**Top row of each group:** difference from the preceding scan. **Bottom row:** each scan’s median pairwise difference to other scans, excluding itself. Thin lines mark session changes. Lower means more similar images. RMS reflects contrast and noise as well as registration; it is not displacement or proof of anatomical accuracy. These images also contributed to the refinement fit.

Identical anatomical masks, per-scan median-intensity normalization, and matched display limits are used before and after refinement. Runs lacking required saved transforms were excluded: 14 for Individual 4 and 16 for Individual 5. Comparisons concern the available acquisitions, not a complete acquisition census. Repeated scans are not independent participants.

## Individual 1

![Matched anatomy-only and within-subject-refined references for Individual 1](assets/coreg-stability/individual_1_matched_alignment.gif)

## Individual 2

![Matched anatomy-only and within-subject-refined references for Individual 2](assets/coreg-stability/individual_2_matched_alignment.gif)

## Individual 3

![Matched anatomy-only and within-subject-refined references for Individual 3](assets/coreg-stability/individual_3_matched_alignment.gif)

## Individual 4

![Matched anatomy-only and within-subject-refined references for Individual 4](assets/coreg-stability/individual_4_matched_alignment.gif)

## Individual 5

![Matched anatomy-only and within-subject-refined references for Individual 5](assets/coreg-stability/individual_5_matched_alignment.gif)

## Individual 6

![Matched anatomy-only and within-subject-refined references for Individual 6](assets/coreg-stability/individual_6_matched_alignment.gif)

[Numerical summaries](assets/coreg-stability/subject_summaries.csv) · [Figure PDF](assets/coreg-stability/six_subject_stability.pdf) · [Longitudinal alignment](longitudinal-alignment.md)

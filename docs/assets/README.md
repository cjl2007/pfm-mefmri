# Curated documentation assets

See the [gallery](../output-gallery.md) for interpretation. `manifest.json`
records file fingerprints and rendering types, not private participant paths.

The README's RidgePFM, crawling-seed FC, coregistration-over-time, and run-level
grayplot examples use ME01. Older summary reportlets retain their original
densely sampled examples. The TANS illustrations retain their separate SAL
example; the workflow diagram is conceptual.

- The coregistration, SDC flicker, and NSI panels are unchanged reportlets.
- The coregistration-over-time GIF cycles through the ME01 anatomical anchor and
  all 44 corrected scan references in chronological order.
- Alignment stability is replotted from saved metrics in numeric session/run
  order, correcting the original display's lexicographic folder order.
- The ME01 motion grayplot shows session 1/run 1 using the production plotting
  code, including scanner-aligned respiratory and pulse traces. The unchanged
  motion-coupling summary covers the older 26-run example.
- The ME01 seed-FC movie uses the production scene with subject-specific
  surfaces and the same concatenated, 1.7 mm-smoothed input used for RidgePFM.
  Its still is decoded at 12 s. The animated preview begins with a three-second
  callout identifying the small white seed, then samples one frame per source
  second and plays at four frames per second.
- Infomap and RidgePFM images are native Workbench documentation renders of
  saved maps. The view layout is adjusted, not the scientific estimates. The
  Infomap view does not borrow RidgePFM network borders.
- The TANS panels reuse the reviewed pfm-tans SAL examples with their existing
  illustrative-placement and hypothetical-dose disclosures.
- The workflow SVG is a conceptual diagram, not participant data.

No raw MRI, CIFTIs, meshes, or operator/session records are distributed here.
Generic titles are not a guarantee of anonymization. Review participant-sharing
and resource permissions before publication; do not add source data to Git.

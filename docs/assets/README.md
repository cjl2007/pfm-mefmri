# Curated documentation assets

See the [gallery](../output-gallery.md) for interpretation. `manifest.json`
records file fingerprints and rendering types, not private participant paths.

The QC, seed-FC, and individual-network assets use the replacement multi-visit
example: 13 sessions and 26 resting-state runs. The TANS illustrations retain
their separate SAL example; the workflow diagram is also retained.

- The coregistration, SDC flicker, and NSI panels are unchanged reportlets.
- The coregistration-over-time GIF cycles through the anatomical anchor and all
  26 saved run panels in numeric order, using the production fade routine with
  fewer transition frames for a compact animation.
- Alignment stability is replotted from saved metrics in numeric session/run
  order, correcting the original display's lexicographic folder order.
- Motion grayplots show session 1/run 1 using the production plotting code and a
  generic title. The unchanged motion-coupling summary covers all 26 runs.
- The seed-FC movie is the original saved movie; its still is decoded at 12 s.
  The animated preview samples one frame per source second and plays at four
  frames per second, looping through the full sweep.
- Infomap and RidgePFM images are native Workbench documentation renders of
  saved maps. The view layout is adjusted, not the scientific estimates. The
  Infomap view does not borrow RidgePFM network borders.
- The TANS panels reuse the reviewed pfm-tans SAL examples with their existing
  illustrative-placement and hypothetical-dose disclosures.
- The workflow SVG is a conceptual diagram, not participant data.

No raw MRI, CIFTIs, meshes, or operator/session records are distributed here.
Generic titles are not a guarantee of anonymization. Review participant-sharing
and resource permissions before publication; do not add source data to Git.

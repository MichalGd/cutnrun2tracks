# Test matrix

The synthetic suite covers config safety, exact control resolution, SE/SEACR
rejection, cross-antibody cohort isolation, consensus biological support,
spike-in formulas/thresholds, checkpoint invalidation/adoption, the indexed
input requirement for canonical BAM filtering, and filtering temporary-cleanup
scope across all four BAM branches. It also covers strict versus continuation
peak-caller failure policies and consensus construction after a failed sample
is excluded from peak contribution. Before release, add
Linux integration fixtures for:

- CUT&RUN/CUT&Tag PE narrow and broad targets;
- CUT&RUN/CUT&Tag SE MACS3 presets;
- all duplicate/MAPQ BAM branches and pair-safe blacklist orphan removal;
- all five coverage families against hand-calculated counts;
- dm6, E. coli, custom, and disabled spike modes;
- multi-condition pairwise differential contrasts and full-rank failures;
- child-process failure, checkpoint resume, cleanup, and UCSC round-trip;
- small real/public CUT datasets for every assay/layout/caller combination.

# Test matrix

The synthetic suite covers config safety, exact control resolution, SE/SEACR
rejection, cross-antibody cohort isolation, consensus biological support,
spike-in formulas/thresholds, checkpoint invalidation/adoption, the indexed
input requirement for canonical BAM filtering, and filtering temporary-cleanup
scope across all four BAM branches. It also covers strict versus continuation
peak-caller failure policies, consensus construction after a failed sample is
excluded from peak contribution, and cohort-local continuation after zero-count
consensus normalization. Nearest-gene annotation is exercised with consensus
and GTF records deliberately presented in different chromosome orders. Unified
reporting is exercised with retained completed-run tables and a synthetic
MultiQC executable; the fixture verifies custom content, parsed data, report
checksums, and that upstream outputs are unchanged. The feature-annotation
fixture covers promoter, enhancer, exon, intron, gene-end, intergenic, and
unclassified assignments; exact absolute/fractional summaries; and plot
creation. Static interface tests cover the self-contained launcher, structured
logging, resource budget, samplesheet/config ownership boundary, cohort
membership manifest, and epic2 auto-primary selection for broad targets.
Before release, add
Linux integration fixtures for:

- CUT&RUN/CUT&Tag PE narrow and broad targets, including epic2 against a
  version-pinned real executable;
- CUT&RUN/CUT&Tag SE MACS3 presets;
- all duplicate/MAPQ BAM branches and pair-safe blacklist orphan removal;
- all five coverage families against hand-calculated counts;
- dm6, E. coli, custom, and disabled spike modes;
- multi-condition pairwise differential contrasts and full-rank failures;
- child-process failure, checkpoint resume, cleanup, and UCSC round-trip;
- small real/public CUT datasets for every assay/layout/caller combination.

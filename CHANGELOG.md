# Changelog

## 0.2.0 - 2026-08-22

- Added a reusable, headless deepTools metagene module for TSS, TES, and scaled
  gene-body profiles with matching PNG/PDF heatmaps.
- Added deterministic GTF-to-BED12 protein-coding reference preparation,
  blacklist/contig/length/overlap filters, and auditable exclusion tables.
- Added HPA broadly expressed reference-set and Ensembl BioMart mouse-ortholog
  builders.
- Added manifest validation, explicit CPM/spike-in track selection, reporting
  tables, resource bounds, and signature-aware workflow integration.

## 0.1.0 - 2026-08-22

- Initial independent CUT&RUN/CUT&Tag workflow implementation.
- Added strict samplesheet/config validation and target-specific cohort IDs.
- Added PE fragment-aware alignment, filtering, coverage, SEACR, and MACS3 paths.
- Added five independent coverage families and generalized spike-in calibration.
- Added target-only differential analysis plus separately labelled control-aware
  sensitivity analyses.

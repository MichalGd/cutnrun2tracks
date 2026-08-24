# Changelog

## 0.2.3 - 2026-08-24

- Scoped each filtering temporary directory and its cleanup trap to a
  subshell, preventing the trap from firing again after its local `tmp`
  variable had expired when a parallel worker returned.
- Added an executable filtering-stage regression fixture that exercises all
  four BAM branches and fails if temporary cleanup leaks into the worker.
- Documented recovery of fully validated filtering outputs without repeating
  completed preprocessing, alignment, or filtering work.

## 0.2.2 - 2026-08-24

- Fixed canonical-contig filtering by applying region selection directly to
  the indexed duplicate-marked BAM instead of an unindexed temporary BAM.
- Reused validated marked BAMs on filtering recovery, avoiding unnecessary
  repeated Picard duplicate marking.
- Made `--from-stage` validate and adopt unchanged earlier-stage checkpoints
  across a workflow signature update, while still forcing the named and later
  stages.
- Added regression coverage for canonical filtering and checkpoint adoption.

## 0.2.1 - 2026-08-24

- Added validated `THREADS_FASTQC` and `THREADS_TRIMGALORE` settings and passed
  them to FastQC and Trim Galore, matching the proven ATACseq2tracks resource
  model for high-capacity servers.
- Made the shared bounded worker pool work-conserving by refilling a slot when
  any worker finishes instead of waiting for submission order.
- Fixed `set -u` failures caused by referencing variables within the same
  `local` declaration in the workflow dispatcher and preprocessing worker.
- Added regression coverage for preprocessing thread validation and parallel
  pool scheduling.

## Unreleased - documentation

- Expanded the GitHub README with a Mermaid workflow diagram, installation,
  stage, output, validation, recovery, and limitations guidance.
- Added a task-oriented documentation index and an exact pipeline-stage guide.
- Clarified output discovery and checkpoint/cleanup recovery without changing
  workflow behavior or configuration defaults.
- Added a server installation runbook covering isolated environment cloning,
  pinned SEACR, shared host-reference reuse, derived canonical/dm6 coordinate
  files, validation, pilot deployment, and guarded promotion.
- Added an explicit all-user deployment model with immutable runtime
  permissions, an optional `/usr/local/bin` launcher, parent-path auditing, and
  acceptance tests executed as an ordinary server user.

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

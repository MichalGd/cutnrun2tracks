# Changelog

## 0.3.0 - 2026-08-30

- Removed duplicated run-wide assay and row-level blacklist settings. Genome,
  layout, assay, FASTQs, biological metadata, and control mappings now live only
  in the samplesheet; assembly references remain only in config.
- Added `cohort_membership.tsv` so targets, controls, shared-control reuse, and
  derived cohort membership are directly auditable.
- Added per-technical-unit FastQC, NRF/PBC library-complexity metrics,
  cohort-specific replicate Spearman/PCA QC, and optional phantompeak
  cross-correlation. preseq remains non-fatal and descriptive.
- Added optional epic2 broad-domain calling through a separate versioned
  environment, with broad/mixed-only validation, matched controls, explicit
  effective-genome settings, caller status, and normal fault isolation.
- Added annotation of every successful per-sample caller peak set and optional
  primary consensus peak set. Outputs include exclusive promoter/enhancer/exon/
  intron/gene-end/other-regulatory/intergenic/unclassified assignments, a full
  overlap table, counts, fractions, base-pair fractions, nearest genes/TSS
  distances, and horizontal stacked PNG/PDF/SVG plots.
- Added bounded parallel pools for spike-in, normalized-track families,
  differential cohorts, annotation, checkpoint hashing, and final checksums.
  Preflight now writes a jobs x threads resource audit and can fail on CPU
  overcommit.
- Added a repository-provided shared launcher that pins the immutable workflow
  release and main Conda environment, allowing a complete run from one
  `cutnrun2tracks --config ...` command without interactive activation.
- Added raw internal console capture, structured workflow/stage/command event
  tables, stage timing, signals/failures, and exit codes while remaining
  compatible with a user-supplied nohup log.
- Fixed strict normalized-track parallel failure handling so the stage writes a
  deterministic FAILED status before returning nonzero.

## 0.2.8 - 2026-08-25

- Treat MultiQC colour-conversion diagnostics as non-fatal when MultiQC exits
  successfully and creates the required report and data directory.
- Retain strict failure handling for non-zero MultiQC exits, module crashes,
  validation errors, and missing report artifacts.
- Add a reporting regression test covering the non-fatal colour diagnostic seen
  in completed human and mouse CUT&Tag report recovery.

## 0.2.7 - 2026-08-25

- Added a final unified MultiQC report over the complete retained workflow
  output, matching the reporting role used by ATACseq2tracks while preserving
  CUT-specific sample, cohort, peak, normalization, and differential semantics.
- Added deterministic MultiQC custom-content tables for retained observations,
  alignment/filtering counts, FRiP, peak calls, consensus sets, normalized-track
  families, spike-in calibration, differential status, and comparison summaries.
- Added selected target/control fingerprints, descriptive TSS profiles, and
  differential PCA/dispersion images as MultiQC custom content while excluding
  the fragile native deepTools parser; authoritative source files are unchanged.
- Made the report stage validate its lightweight HTML/TSV outputs and its final
  MultiQC HTML before cleanup can proceed.
- Added `utilities/regenerate_reports.sh` to recover the lightweight and unified
  reports from completed analyses without rerunning upstream stages. A separate
  report checksum inventory avoids rewriting historical workflow checksums.
- Expanded the lightweight report to include normalization, differential, and
  preseq warnings, and added end-to-end reporting regression coverage.

## 0.2.6 - 2026-08-25

- Fixed nearest-gene annotation failures caused by inconsistent chromosome
  ordering between consensus and GTF-derived BED inputs.
- Genome-sorted both inputs with the configured chromosome-sizes file and
  passed that same order to `bedtools closest -g`, making annotation independent
  of shell locale and suitable for hg38, mm39, and compatible custom genomes.
- Added an executable annotation-order regression test.

## 0.2.5 - 2026-08-25

- Made consensus-normalization failures cohort-local when
  `REQUIRE_ALL_ENABLED_TRACKS=false`, so a zero-count or otherwise
  non-normalizable cohort is recorded and skipped without terminating
  unaffected cohorts, QC, metagene, annotation, or reporting.
- Kept normalization scientifically strict: samples with zero consensus counts
  are not silently removed and no scaling factors are fabricated. Strict
  fail-fast behavior remains available with `REQUIRE_ALL_ENABLED_TRACKS=true`.
- Added per-sample `consensus_count_sums.tsv`, per-family status tables, and
  dedicated factor-calculation logs to make normalization exclusions auditable.
- Made differential analysis recognize an upstream normalization skip as an
  expected cohort skip while continuing to treat unexpectedly missing count
  tables as failures.
- Added continuation/strict-mode regression coverage with one failing and one
  successful cohort.

## 0.2.4 - 2026-08-24

- Added configurable `PEAKCALL_FAILURE_POLICY=continue|fail`; continuation mode
  records per-caller `SUCCESS`, `EMPTY`, or `ERROR` results without allowing one
  failed sample or auxiliary caller to terminate all other samples.
- Made consensus construction exclude only failed/empty primary peak
  contributions, proceed when enough successful biological samples remain, and
  record every exclusion and reason in machine-readable TSV files.
- Preserved all BAM, coverage, QC, metagene, and downstream consensus-counting
  eligibility for a peak-call-excluded sample; peak failure is reported and is
  not treated as an automatic biological QC exclusion.
- Added completed-with-warnings reporting and regression coverage for both
  continuation and strict peak-calling policies.

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

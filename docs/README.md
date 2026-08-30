# cutnrun2tracks documentation

This index groups the documentation by user task. Start with the root
[README](../README.md) for installation, quick start, and the principal output
tree. The pages below describe version 0.3.1 unless stated otherwise.

## Prepare and run an analysis

| Topic | Document |
|---|---|
| Launch an analysis on the established shared server | [One-page shared-server quick start](19_server_quick_start.md) |
| Samplesheet rows, technical/biological replicates, controls, and validation | [Inputs and configuration](01_inputs_and_configuration.md) |
| Plan technical/biological replicates, matched controls, cohorts, blocks, and contrasts | [Replicates and experimental design](18_replicates_and_experimental_design.md) |
| Complete configuration keys and example samplesheets | [`config/README.md`](../config/README.md) and [`config/config.conf.template`](../config/config.conf.template) |
| Stage order, optional stages, checkpoints, and declared outputs | [Pipeline stages](07_pipeline_stages.md) |
| Shared server installation, launchers, environments, and exact deployed references | [Server installation](08_server_installation.md) |
| Build, validate, and record every host, spike-in, annotation, and metagene reference | [Reference preparation and provenance](17_reference_preparation_and_provenance.md) |
| Move an older project to the 0.3 metadata contract | [Migration from 0.2](14_migration_from_0.2.md) |

## Methods and interpretation

| Topic | Document |
|---|---|
| Overall analytical flow and signal-unit definitions | [Methods overview](02_methods.md) |
| SEACR, MACS3, epic2, parameters, controls, failure handling, and consensus construction | [Peak calling and consensus](09_peak_calling.md) |
| Canonical-contig, mapping-quality, duplicate, mitochondrial, and exact blacklist filtering | [Reference filtering and blacklists](10_references_blacklist_and_filtering.md) |
| FastQC/trimming/alignment QC and every assay-performance metric | [Quality control](11_quality_control.md) |
| CPM, DESeq2-consensus, robust-CPM, and spike-in track families and formulas | [Tracks and spike-in normalization](12_tracks_and_normalization.md) |
| Primary statistical models, controls, blocking, contrasts, and sensitivity analyses | [Differential binding](03_differential_enrichment.md) |
| Consensus, differential-result, and all-peak feature annotation | [Genomic annotation](13_genomic_annotation.md) |
| TSS/TES/gene-body aggregate plots and gene-set references | [Metagene aggregate-signal module](06_metagene.md) |

## Operate and validate the workflow

| Topic | Document |
|---|---|
| Diagnose failures, warnings, skips, resource issues, and safe restart stages | [Troubleshooting](16_troubleshooting.md) |
| Output organization, logs, checkpoint recovery, partial reruns, and cleanup | [Outputs and recovery](04_outputs_and_recovery.md) |
| Scientific limitations and pilot decisions | [Limitations](05_limitations.md) |
| Automated test coverage and outstanding real-data validation | [Test matrix](06_test_matrix.md) |
| Compare this manual with ATACseq2tracks and prioritize missing pages | [Documentation coverage comparison](15_documentation_comparison.md) |

## Project information

- [Release history](../CHANGELOG.md)
- [Contribution requirements](../CONTRIBUTING.md)
- [Metagene shared-module interface](../common/metagene/README.md)

Release history and retired behavior are kept outside the current-use path in
`CHANGELOG.md` and the migration page. Method pages describe the installed
0.3.1 behavior unless they explicitly identify a limitation or historical
compatibility note.

Executable behavior is defined by `cutnrun2tracks.sh`, `scripts/`, and the
validated configuration template. For a completed run, inspect
`00_metadata/resolved_config.tsv`, `sample_manifest.tsv`,
`reference_manifest.tsv`, and `software_versions.tsv` before assuming that a
template default was used.

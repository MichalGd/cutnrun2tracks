# cutnrun2tracks

`cutnrun2tracks` is a samplesheet-driven workflow for paired-end and single-end
CUT&RUN and CUT&Tag data. One command validates the experiment, processes all
libraries, calls and consolidates peaks, creates normalized browser tracks,
runs quality control and differential enrichment, annotates peaks, and writes
final HTML reports.

Current release: **0.3.1**. The tagged Linux release and its main and optional
epic2 environments have passed the repository test suite. Scientific
interpretation still requires experiment-specific QC and appropriate biological
replication.

## What the workflow provides

- technical-replicate FASTQ merging while biological replicates remain
  independent;
- exact replicate- and context-matched IgG, input, or mock controls;
- FastQC, Trim Galore, Bowtie2 alignment, duplicate marking, canonical-contig,
  blacklist, MAPQ, and duplicate-policy filtering;
- four inspectable filtering branches and a samplesheet-selected analysis BAM;
- CPM, DESeq2-consensus, robust-CPM, and optional spike-in-calibrated tracks;
- matched-control MACS3 and SEACR peak calling plus optional epic2 broad-domain
  calling, with per-sample failure isolation;
- caller-, peak-class-, target-, and antibody-specific consensus peaks based on
  biological-sample support;
- fragment/read counts, FRiP, fragment length, NRF/PBC, preseq, fingerprints,
  replicate correlation/PCA, descriptive TSS QC, and optional cross-correlation;
- optional TSS, TES, and scaled-gene-body aggregate-signal plots;
- primary target-only DESeq2Enrichment and DiffBind analyses plus clearly
  separated control-aware sensitivity models;
- annotation of every successful per-sample peak set and primary consensus into
  promoter, enhancer, exon, intron, gene-end, other-regulatory, intergenic, or
  unclassified categories, including all-overlap tables and horizontal stacked
  plots; and
- structured logs, content-aware checkpoints, browser assets, a unified
  MultiQC report, a lightweight HTML report, provenance, and final checksums.

The workflow is an independent CUT&RUN/CUT&Tag implementation informed by the
engineering patterns in
[ATACseq2tracks](https://github.com/MichalGd/ATACseq2tracks). It does not call or
modify an ATACseq2tracks installation.

## Run an analysis on the shared server

The system launcher pins the immutable workflow release and Conda environment;
users do not activate Conda or export tool paths manually.

Prepare one configuration file and the samplesheet referenced by its
`SAMPLESHEET` key, then run:

```bash
cutnrun2tracks --config /absolute/path/to/config.conf
```

For the installed release explicitly:

```bash
cutnrun2tracks-0.3.1 --config /absolute/path/to/config.conf
```

Before processing reads, validate metadata and resources:

```bash
cutnrun2tracks --config /absolute/path/to/config.conf --plan
cutnrun2tracks --config /absolute/path/to/config.conf --preflight-only
```

The normal full run needs no additional arguments. A simple external console
log can be retained with:

```bash
PROJECT=/absolute/path/to/project
OUTPUT="$PROJECT/cutntag2tracks"
CONFIG="$PROJECT/config/config.conf"

mkdir -p "$OUTPUT"
nohup cutnrun2tracks --config "$CONFIG" \
    > "$OUTPUT/cutnrun2tracks.nohup.log" 2>&1 &
echo $! > "$OUTPUT/cutnrun2tracks.pid"
```

The workflow also writes its own raw console log and structured event tables,
so the `nohup` file is a convenient second capture rather than the only log.

## Inputs: one source of truth per field

Start from:

```text
config/config.conf.template
config/samplesheet_template.csv
config/examples/
```

The division of responsibility is intentional:

| File | Contains |
|---|---|
| `config.conf` | samplesheet and output paths, genome-reference paths, workflow policies, thresholds, optional modules, and resource limits |
| `samplesheet.csv` | genome, PE/SE layout, assay profile, FASTQ paths, technical and biological replicate identity, factor/antibody, condition, target/control role, and control mapping |

Do not duplicate row-level genome, layout, assay, FASTQ, or biological metadata
in config. One samplesheet row is one sequencing unit. Rows with the same
`sample_id` and biological `replicate` but different `tech_replicate` values are
merged before trimming; different biological replicates are never merged.

The validator resolves target-control mappings and cohorts before any expensive
work. Inspect:

```text
<OUTPUT_DIR>/00_metadata/sample_manifest.tsv
<OUTPUT_DIR>/00_metadata/cohort_membership.tsv
<OUTPUT_DIR>/00_metadata/reference_manifest.tsv
<OUTPUT_DIR>/00_metadata/resource_budget.tsv
```

See [Inputs and configuration](docs/01_inputs_and_configuration.md) for the
complete schema and worked examples.

## Analytical outline

```text
FASTQ + metadata
  -> preflight
  -> FastQC, merge, trim
  -> host alignment and optional competitive spike-in alignment
  -> duplicate marking and filtering branches
  -> CPM tracks and peak calling
  -> biological-support consensus
  -> spike-in and DESeq2-derived track families
  -> metagene and assay QC
  -> differential enrichment
  -> gene/cCRE/genomic-feature annotation and browser assets
  -> unified reports, guarded cleanup, and final checksums
```

The exact stage order is documented in
[Pipeline stages](docs/07_pipeline_stages.md). Long-running stages use bounded
parallel pools, and preflight checks jobs x threads against `TOTAL_CPU_BUDGET`.

## Controls, cohorts, and normalization

Matched controls have three separate roles:

1. background for peak calling;
2. target-control QC; and
3. optional, explicitly labelled sensitivity analyses.

The primary differential model counts raw target fragments or reads over a
condition-unbiased cohort consensus. Controls are not biological replicates,
are not subtracted from primary DESeq2 input, and do not determine its size
factors. Spike-in calibration is a separate operation.

A cohort includes genome, assay profile, factor, antibody ID, layout, target
class, duplicate policy, primary caller, and primary peak class. Different
targets or antibody IDs are never normalized together. Consensus support is
counted across biological sample keys after technical merging; the default
minimum of two is a pragmatic reproducibility rule, not formal IDR.

Read:

- [Differential binding](docs/03_differential_enrichment.md) for models,
  controls, contrasts, blocking, and sensitivity analyses;
- [Tracks and spike-in normalization](docs/12_tracks_and_normalization.md) for
  every track family and formula; and
- [Peak calling and consensus](docs/09_peak_calling.md) for caller-specific
  parameters and fault handling.

## Principal outputs

```text
<OUTPUT_DIR>/
|-- 00_metadata/     resolved metadata, provenance, resource audit, event logs
|-- 01_fastq_qc/     raw and trimmed FastQC/MultiQC
|-- 02_trimmed_fastq/
|-- 03_alignment/    marked, filtered, analysis, and optional spike-in BAMs
|-- 04_tracks/       CPM, DESeq2, robust-CPM, and spike-in tracks
|-- 05_peaks/        per-sample caller output and cohort consensuses
|-- 06_qc/           assay QC, replicate QC, and optional metagene results
|-- 07_annotation/   nearest-gene, cCRE, all-peak tables, and composition plots
|-- 08_differential/ primary and sensitivity analyses
|-- 09_browser/      UCSC track definitions and IGV session
|-- 10_reports/      unified MultiQC and lightweight HTML reports
|-- logs/            workflow and tool logs
`-- .checkpoints/    signature-and-output stage checkpoints
```

The main report is:

```text
10_reports/cutnrun2tracks_multiqc_report.html
```

See [Outputs and recovery](docs/04_outputs_and_recovery.md) for complete file
contracts and cleanup behavior, and [Genomic annotation](docs/13_genomic_annotation.md)
for annotation sources, columns, classification rules, and limitations.

## Resume and troubleshoot a run

Every completed stage has a content-aware JSON checkpoint. To validate all
earlier stages and rerun a named stage plus everything after it:

```bash
cutnrun2tracks \
    --config /absolute/path/to/config.conf \
    --from-stage qc
```

To stop cleanly after a stage:

```bash
cutnrun2tracks \
    --config /absolute/path/to/config.conf \
    --stop-after consensus
```

The internal logs are:

```text
logs/cutnrun2tracks.console.log
00_metadata/workflow_events.tsv
00_metadata/command_events.tsv
00_metadata/stage_timing.tsv
```

Do not copy or edit checkpoint JSON to bypass validation. If a code or metadata
change can affect an earlier output, rerun from that earlier stage. Set
`ENABLE_AUTOMATIC_CLEANUP=false` before a run when intermediate FASTQs/BAMs are
needed for planned troubleshooting.

Reports can be rebuilt from retained completed-run outputs without rerunning
upstream analysis:

```bash
bash utilities/regenerate_reports.sh --output-dir /absolute/path/to/results
```

## Installation and validation

Administrators should use the versioned release, launcher, main environment,
and optional epic2 sidecar procedure in
[Server installation](docs/08_server_installation.md). A local installation can
be created from `environment.yml`; SEACR and epic2 remain separately pinned as
described in that runbook.

Repository validation in the locked Linux environment:

```bash
bash tests/check_bash_syntax.sh
bash tests/run_tests.sh
```

Passing synthetic tests establishes implementation consistency, not biological
validity for a new protocol or antibody. Review the final QC report, excluded
peak samples, consensus status, normalization-family status, and differential
model diagnostics before interpreting tracks or peaks.

## Documentation

| Question | Page |
|---|---|
| How do I launch on the shared server? | [One-page server quick start](docs/19_server_quick_start.md) |
| How do I prepare metadata and controls? | [Inputs and configuration](docs/01_inputs_and_configuration.md) |
| How should replicates, controls, batches, and donors be designed? | [Replicates and experimental design](docs/18_replicates_and_experimental_design.md) |
| What exactly does each analysis step do? | [Methods](docs/02_methods.md) and [Pipeline stages](docs/07_pipeline_stages.md) |
| How are peaks called? | [Peak calling and consensus](docs/09_peak_calling.md) |
| Which QC metrics are produced? | [Quality control](docs/11_quality_control.md) |
| How are tracks and spike-in factors calculated? | [Tracks and normalization](docs/12_tracks_and_normalization.md) |
| How does differential binding work? | [Differential binding](docs/03_differential_enrichment.md) |
| Which annotation sources and rules are used? | [Genomic annotation](docs/13_genomic_annotation.md) |
| How are all reference families built and audited? | [Reference preparation and provenance](docs/17_reference_preparation_and_provenance.md) |
| How do I recover outputs or reports? | [Outputs and recovery](docs/04_outputs_and_recovery.md) |
| How do I diagnose a failed or skipped stage? | [Troubleshooting](docs/16_troubleshooting.md) |
| How do I migrate a 0.2 project? | [Migration from 0.2](docs/14_migration_from_0.2.md) |
| What documentation is still missing? | [Documentation comparison](docs/15_documentation_comparison.md) |

The complete task-oriented index is [docs/README.md](docs/README.md).
Release history belongs in [CHANGELOG.md](CHANGELOG.md), not in the quick-start
path. Scientific and implementation limits are collected in
[docs/05_limitations.md](docs/05_limitations.md).

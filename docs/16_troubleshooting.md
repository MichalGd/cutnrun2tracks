# Troubleshooting

[Documentation index](README.md) | [Outputs and recovery](04_outputs_and_recovery.md)

This guide applies to `cutnrun2tracks` 0.3.1. Diagnose the failed or skipped
unit before rerunning anything. A workflow-wide failure, a deliberately skipped
cohort, an empty biological peak call, and a non-fatal QC warning have different
meanings.

## Five-minute triage

Set the project paths once:

```bash
PROJECT=/absolute/path/to/project
OUTPUT="$PROJECT/cutntag2tracks"
CONFIG="$PROJECT/config/config.conf"
```

### 1. Is the workflow still running?

```bash
if [[ -s "$OUTPUT/cutnrun2tracks.pid" ]]; then
    ps -fp "$(cat "$OUTPUT/cutnrun2tracks.pid")" || true
fi

pgrep -af 'cutnrun2tracks|fastqc|trim_galore|cutadapt|bowtie2|samtools|macs3|SEACR|epic2|bamCoverage|computeMatrix|preseq' || true
```

A quiet parent process can be waiting for active workers. Check child tools and
logs before concluding that the run is stuck.

### 2. Which stage last started or completed?

```bash
grep '^=== \[' "$OUTPUT/logs/cutnrun2tracks.console.log" | tail -n 25

if [[ -s "$OUTPUT/00_metadata/workflow_events.tsv" ]]; then
    tail -n 25 "$OUTPUT/00_metadata/workflow_events.tsv"
fi

if [[ -s "$OUTPUT/00_metadata/stage_timing.tsv" ]]; then
    column -t -s $'\t' "$OUTPUT/00_metadata/stage_timing.tsv"
fi
```

The last `START` without `COMPLETE`, `SKIPPED`, or `REUSED` identifies the
stage to inspect.

### 3. Read both console captures

```bash
tail -n 100 "$OUTPUT/logs/cutnrun2tracks.console.log"
tail -n 100 "$OUTPUT/cutnrun2tracks.nohup.log" 2>/dev/null || true
```

The internal console log is authoritative. A user-named `nohup` log is an
additional raw capture. When using `tail -f`, do not append an accidental `~`
to the filename.

### 4. Find explicit failures and skips

```bash
find "$OUTPUT" -type f \
    \( -name 'FAILED.json' -o -name 'SKIPPED.json' -o -name '*status.tsv' \) \
    -print | sort

grep -RniE '(^|[^A-Z])(ERROR|FAILED|WARNING|SKIPPED)([^A-Z]|$)' \
    "$OUTPUT/logs" "$OUTPUT/00_metadata" 2>/dev/null | tail -n 100
```

Read a JSON/status file in its cohort or sample context. Do not treat every
`SKIPPED` as an error.

### 5. Check checkpoints

```bash
find "$OUTPUT/.checkpoints" -maxdepth 1 -type f -name '*.json' -printf '%f\n' | sort
```

A checkpoint exists only after a stage and its declared outputs have completed
validation. Its presence does not authorize manually copying it to another run.

## Status vocabulary

| Status | Meaning | Normal response |
|---|---|---|
| `SUCCESS` / `COMPLETE` | required work and output validation completed | continue |
| `REUSED` | earlier outputs passed content-aware validation during `--from-stage` | continue |
| `SKIPPED` | module/cohort was disabled or ineligible by a documented rule | review reason; rerun only if inputs/policy should change |
| `EMPTY` | caller ran but produced no valid peaks | review sample QC; not automatically a software error |
| `ERROR` | a per-sample peak caller failed | inspect caller log; continuation may protect other samples |
| `WARNING` | non-fatal diagnostic, such as preseq/cross-correlation failure | review but do not assume workflow failure |
| `FAILED` | required command, validation, strict policy, or stage failed | correct cause, then resume from earliest affected stage |

## Input and preflight problems

### Unknown or obsolete config key

`CONFIG ERROR: unknown key` usually means a config from an older release was
reused. Start from the installed 0.3.1 template and copy project values; do not
delete validation rules or source the config as shell code. See
[Migration from 0.2](14_migration_from_0.2.md).

### Windows CRLF or encoding correction

The input sanitizer may report `Corrected ... (CRLF)` and create a timestamped
`*.windows-artifact-backup.*`. This is expected recovery, not a failure. Review
the normalized file and keep the backup until validation succeeds.

### Samplesheet header mismatch

The header must exactly match `config/samplesheet_template.csv`, including
order. Do not add a row-level blacklist or remove empty spike-in columns.

### Target cannot resolve its control

Typical messages include `expected exactly one matched control`, `incompatible
in ...`, or a control assigned to multiple targets. Check:

- `control_id` equals the control row's `sample_id`, case-sensitively;
- target and ordinary control have the same biological `replicate`;
- genome, assay profile, layout, treatment, cell type, batch, and spike metadata
  match; and
- condition also matches unless the deliberately shared-control policy is used.

Run `--plan` and inspect `cohort_membership.tsv`; do not repair control mapping
by relabelling biological replicates as technical replicates.

### Missing reference or index

Preflight reports the exact missing key/path. Verify readability, non-empty
content, Bowtie2 prefix components, matching assembly, and chromosome naming.
Use [Reference preparation and provenance](17_reference_preparation_and_provenance.md),
then rerun `--preflight-only`.

### Resource-budget failure

When `RESOURCE_CHECK_MODE=fail`, preflight stops if jobs x threads exceeds
`TOTAL_CPU_BUDGET`. Inspect:

```bash
column -t -s $'\t' "$OUTPUT/00_metadata/resource_budget.tsv"
```

Reduce the relevant `*_PARALLEL_JOBS`, its per-tool thread count, or both. Do
not increase the CPU budget above the cores actually available to the job.

### epic2 smoke check or `pkg_resources` failure

The optional epic2 0.0.52 sidecar requires the version-pinned launcher and
`setuptools=80.9.0`. Confirm the installed workflow and sidecar:

```bash
cutnrun2tracks --version
/usr/local/bin/epic2-0.0.52 --help >/dev/null
/opt/miniconda/bin/conda list \
  --prefix /opt/miniconda/envs/cutnrun2tracks-epic2-0.3.0 |
  grep -E '^(epic2|setuptools)[[:space:]]'
```

If epic2 is not required, remove it from `PEAK_CALLERS`; do not point at an
unvalidated environment.

## Stage-specific diagnostics

### Preprocessing appears under-parallelized

The worker pool processes libraries concurrently, while FastQC and Trim Galore
use their own per-worker thread counts. Compare live processes with
`QC_SAMPLE_PARALLEL_JOBS`, `THREADS_FASTQC`, and `THREADS_TRIMGALORE`. CPU use
can be limited by gzip decompression, storage throughput, short files, or a
phase in which some workers have already finished. The resource-budget table
shows the configured ceiling, not a guarantee of constant 100% CPU occupancy.

Check progress by counting completed merged/trimmed pairs and reading the most
recent per-sample logs. Do not launch a second preprocessing run into the same
output directory.

### Alignment or BAM integrity failure

For retained BAMs:

```bash
SAMTOOLS=/opt/miniconda/envs/cutnrun2tracks-0.3.0/bin/samtools
find "$OUTPUT/03_alignment" -type f -name '*.bam' -print0 |
    xargs -0 -r "$SAMTOOLS" quickcheck
```

An empty/truncated BAM, missing index, or failed `idxstats` means the stage is
not reusable. Investigate disk space, the Bowtie2 log, input FASTQ integrity,
and the exact command event before rerunning alignment.

### Filtering removed every signal unit

Inspect `03_alignment/metrics/*.filter_counts.tsv`, alignment rates, MAPQ,
canonical contigs, blacklist assembly, proper-pair filtering, and duplicate
policy. Do not set `ALLOW_EMPTY_FILTERED_BAM=true` merely to make the workflow
continue; first establish whether the sample is biologically empty or the
reference/filter contract is wrong.

### Peak caller failed for one sample

Read:

```text
05_peaks/per_sample/<sample>/caller_status.tsv
05_peaks/per_sample/peakcall_status.tsv
logs/peakcalling/<sample>.*.log
```

With `PEAKCALL_FAILURE_POLICY=continue`, one caller/sample failure is recorded
and unaffected samples continue. `fail` makes any caller error fatal. SEACR can
fail or produce no peaks for extremely sparse target/control bedGraphs; MACS3
or epic2 output for the same sample is independent. An `awk: division by zero`
inside SEACR is a caller/sample failure, not evidence that all other samples
should be discarded.

### Consensus is skipped

Read `05_peaks/consensus/consensus_summary.tsv` and the cohort's
`excluded_peak_samples.tsv`. Common reasons are:

- too few successful primary peak samples after exclusions;
- fewer than `CONSENSUS_MIN_BIOLOGICAL_SAMPLES` independent samples; or
- no intervals meeting the support threshold.

Technical FASTQ units never count as independent support. Lowering support or
allowing a single-sample consensus changes the scientific definition and
requires explicit justification.

### Normalized-track family is skipped

Read:

```text
04_tracks/normalized_track_family_status.tsv
04_tracks/<family>/<cohort>/tables/consensus_count_sums.tsv
logs/normalized_tracks/<cohort>.*.factors.log
```

Zero counts for one or more samples, unavailable consensus, or failed factor
estimation can invalidate only that cohort/family. With
`REQUIRE_ALL_ENABLED_TRACKS=false`, other cohorts and later stages continue.
Strict mode deliberately fails the stage.

### `preseq` failed

`preseq` is descriptive and non-fatal. Messages such as `max count before zero
is less than min required count` occur when the retained library has too few
duplicate-count levels to fit an extrapolation. Review NRF/PBC and depth; the
absence of a preseq curve does not invalidate otherwise intact BAMs.

### Metagene stage cannot resolve tasks

Verify `METAGENE_GENE_SET_MANIFEST`, every BED12 path/checksum, selected gene
set IDs, genome labels, track-family availability, blacklists, and `pyBigWig`.
If spike-in tracks were requested but unavailable, either repair spike-in
normalization or explicitly permit a CPM fallback; do not relabel CPM as
spike-normalized.

### Differential analysis is skipped

Read the relevant `SKIPPED.json`, comparison summary, and
`04_tracks/deseq2_consensus/<cohort>/tables/raw_counts.tsv.gz`. Common valid
reasons are fewer than two replicated conditions, incomplete controls for a
sensitivity model, a block requested in the DiffBind wrapper, or an unavailable
consensus. The block-aware direct DESeq2 path remains the primary analysis.

Design-rank errors require fixing confounding or metadata, not deleting a term
after inspecting results. See
[Replicates and experimental design](18_replicates_and_experimental_design.md).

### Annotation fails or looks implausible

Confirm GTF, cCRE, chromosome sizes, and peaks share an assembly and chromosome
vocabulary. Inspect `peak_annotation_status.tsv`, invalid/unclassified counts,
and the run reference manifest. Large `unclassified` fractions usually indicate
contigs absent from chromosome sizes; unexpected category fractions can reflect
the cCRE release, peak widths, or configured precedence. See
[Genomic annotation](13_genomic_annotation.md).

### MultiQC reports a colour diagnostic

An `mqc_colour` conversion message is non-fatal when MultiQC exits zero and the
HTML plus data directory exist. Version 0.2.8 and later tolerate this diagnostic
but still fail on a non-zero MultiQC exit, parser/module crash, missing report,
or validation error. Check `10_reports/multiqc_status.tsv`.

Rebuild reports without upstream analysis:

```bash
bash /opt/bioinformatics/workflows/cutnrun2tracks/current/utilities/regenerate_reports.sh \
    --output-dir "$OUTPUT"
```

### Disk space or filesystem write failure

Stop new work, check `df -h` and quota, and preserve logs. An incomplete output
must not be accepted by manually writing a checkpoint. After space is restored,
rerun from the earliest stage that may have a partial file. Automatic cleanup
runs only after reporting succeeds; disable it in advance when intermediates
are needed for diagnosis.

## Safe stopping and resuming

Verify the PID belongs to this exact output/config before stopping it:

```bash
PIDFILE="$OUTPUT/cutnrun2tracks.pid"
RUN_PID="$(cat "$PIDFILE")"
ps -fp "$RUN_PID"
kill "$RUN_PID"
```

Do not use broad `pkill` patterns on a multi-user server.

After correcting the cause, choose the earliest affected stage:

| Changed or failed item | Earliest normal restart |
|---|---|
| samplesheet, FASTQs, trimming policy | `preprocess` (often a new output directory for metadata changes) |
| host index/alignment settings | `alignment` |
| blacklist, canonical contigs, MAPQ, duplicate policy | `filtering` |
| coverage settings only | `cpm` |
| caller parameters or peak policy | `peakcalling` |
| consensus support policy | `consensus` |
| spike-in factors/references | `spikein` |
| normalized-track settings | `normalized_tracks` |
| gene-set manifest/metagene settings | `metagene` |
| QC settings | `qc` |
| model, blocks, contrasts, thresholds | `differential` |
| GTF/cCRE/annotation settings | `annotation` |
| reporting only | report-regeneration utility or `report` |

Resume with validation of every earlier checkpoint:

```bash
cutnrun2tracks \
    --config "$CONFIG" \
    --from-stage <stage>
```

Changing a reference or samplesheet can alter earlier scientific outputs; do
not select a later stage merely to save time. Never edit checkpoint hashes or
copy checkpoints between projects.

## Information to retain when requesting help

- workflow version and commit (`cutnrun2tracks --version`);
- config and sanitized samplesheet, with sensitive paths redacted if needed;
- `run_manifest.tsv`, `sample_manifest.tsv`, `cohort_manifest.tsv`,
  `cohort_membership.tsv`, `reference_manifest.tsv`, and `resource_budget.tsv`;
- internal console log, external `nohup` log, workflow/command events, and the
  failing tool's per-sample log;
- failed/skipped status files and the last valid checkpoint;
- disk/quota state and relevant process listing; and
- the exact rerun command.

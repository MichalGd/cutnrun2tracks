# Outputs and recovery

Results use numbered directories from `00_metadata` through `10_reports`.
Sample files use `<output_prefix>.bioR<replicate>`. Cohort IDs contain a readable
prefix and an eight-character SHA-256 suffix; the complete key is retained in
`cohort_manifest.tsv`.

Per-caller results are recorded in
`05_peaks/per_sample/<sample_key>/caller_status.tsv`; the stable run-wide table
is `05_peaks/per_sample/peakcall_status.tsv`. A continuation-mode exclusion from
consensus is recorded under each cohort as `excluded_peak_samples.tsv`. These
records distinguish a biological zero-peak result (`EMPTY`) from a caller
execution failure (`ERROR`) and must be reviewed with the ordinary QC outputs.

## Output map

```text
<OUTPUT_DIR>/
├── 00_metadata/                    run, sample, cohort and reference manifests
├── 01_fastq_qc/                    raw/trimmed FastQC and MultiQC
├── 02_trimmed_fastq/               merged/trimmed FASTQs; removed by default
├── 03_alignment/                   sorted, marked, filtered and analysis BAMs
├── 04_tracks/                      CPM, DESeq2-derived and spike-in tracks
├── 05_peaks/                       per-sample and consensus peaks
├── 06_qc/                          alignment, fragment, FRiP, control and metagene QC
├── 07_annotation/                  consensus and differential annotations
├── 08_differential/                primary and control-aware sensitivity results
├── 09_browser/                     UCSC and IGV assets
├── 10_reports/                     pipeline HTML report
├── logs/                           stage/tool logs
└── .checkpoints/                   signature-and-output JSON checkpoints
```

The root README gives a more detailed
[principal output tree](../README.md#principal-outputs), and
[Pipeline stages](07_pipeline_stages.md) maps every stage to its declared
checkpoint output.

## Checkpoints

Each stage checkpoint is JSON containing the complete run signature and hashes
of declared outputs. The signature covers the sanitized CSV, resolved config,
workflow version, scripts, reusable `common/` modules, and reference manifest.
A changed or missing output invalidates the checkpoint.

Run `--plan` to validate the metadata model and inspect
`00_metadata/planned_stages.tsv` without checking files/tools. Run
`--preflight-only` for the complete input, tool, and reference audit.

Reuse validated outputs before a stage, then force that stage and all later
stages:

```bash
bash cutnrun2tracks.sh --config /path/to/config.conf --from-stage qc
```

This is an explicit recovery override. Before the requested stage, each existing
checkpoint is validated against its stored output sizes and SHA-256 hashes and
its signature adoption is recorded in the checkpoint JSON. If any earlier
checkpoint or output is missing or changed, recovery stops instead of silently
rerunning or accepting it. The named stage and all later stages always rerun.
Select a starting stage at or before the earliest output that could be affected
by the configuration or code change.

If a stage command completes every output but fails during terminal bookkeeping
before its checkpoint is written, preserve the outputs and validate them with
the relevant native tools before recovery. A replacement checkpoint may be
written with `scripts/checkpoint.py write` only after expected file counts,
indexes, links, and content integrity have all been confirmed and the incident
has been recorded. Deploy the code correction, then use `--from-stage` on the
next stage so the replacement checkpoint is explicitly adopted into the new
workflow signature. Never use this procedure for partial or merely
size-checked outputs.

Stop after a stage without editing the workflow:

```bash
bash cutnrun2tracks.sh --config /path/to/config.conf --stop-after consensus
```

## Stable reporting interfaces

Metagene outputs are retained under `06_qc/metagene`. Each sample/gene-set/mode
directory contains the deepTools matrix, exported profile values, sorted BED,
PNG/PDF profile, PNG/PDF heatmap, and task metadata. The run-wide
`artifacts.tsv` is the stable reporting interface.

The final report is `10_reports/pipeline_report.html`. Machine-readable TSVs in
`00_metadata`, `05_peaks/consensus`, `06_qc`, and `08_differential` should be
preferred over parsing HTML or filenames. Browser definitions are written under
`09_browser/`.

## Cleanup and retention

Cleanup begins only after the final report exists, targets only explicit
children of the resolved output root, and records every deletion in
`00_metadata/cleanup_manifest.tsv`.

Set `ENABLE_AUTOMATIC_CLEANUP=false` before starting a run that must retain all
intermediates. The `KEEP_*` switches control individual intermediate families.
Final tracks, peaks, QC, differential results, annotations, browser assets,
reports, metadata, logs, and checksums are not cleanup targets under the default
policy.

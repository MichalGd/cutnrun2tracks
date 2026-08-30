# Shared-server quick start

[Documentation index](README.md) | [Full server installation](08_server_installation.md)

This page is for users of the existing shared `biolserv` deployment. It assumes
the administrator has installed and validated `cutnrun2tracks` 0.3.1, the main
environment, optional sidecars, and shared references. No interactive Conda
activation or PATH export is required.

## 1. Confirm the installed launcher

```bash
cutnrun2tracks --version
cutnrun2tracks --help | head -n 30
```

Expected version:

```text
cutnrun2tracks 0.3.1
```

If another version is current, use its matching documentation/template or the
explicit versioned launcher; do not combine a new executable with an old config
without migration.

## 2. Create a project configuration directory

```bash
PROJECT=/home/$USER/Analysis/my_cut_project
CONFIG_DIR="$PROJECT/config"
OUTPUT="$PROJECT/cutntag2tracks"
RELEASE=/opt/bioinformatics/workflows/cutnrun2tracks/current

mkdir -p "$CONFIG_DIR" "$OUTPUT"
cp "$RELEASE/config/config.conf.template" "$CONFIG_DIR/config.conf"
cp "$RELEASE/config/samplesheet_template.csv" "$CONFIG_DIR/samplesheet.csv"
```

Use a new output directory for a new scientific run. Do not place project
configuration or results inside the root-owned release tree.

## 3. Prepare the samplesheet

One row is one sequencing unit. Enter absolute FASTQ paths and preserve the
exact header. Important rules:

- PE rows need both FASTQs; SE rows require an empty `fastq_2`;
- technical lanes share `sample_id + replicate` and use distinct
  `tech_replicate` values;
- biological replicates use distinct `replicate` values;
- targets use `control_id` to name their matched control `sample_id`;
- targets use `control_type=none`; controls use
  `target_class=control` and `control_type=igg|input|mock`; and
- leave spike fields empty when `SPIKEIN_MODE=none`.

For a starting example instead of the empty template:

```bash
cp "$RELEASE/config/examples/cuttag_pe.csv" "$CONFIG_DIR/samplesheet.csv"
```

Replace every illustrative `/data/...` path and sample value. See
[Replicates and experimental design](18_replicates_and_experimental_design.md)
before encoding shared controls, technical replicates, blocking, or spike-in.

## 4. Edit `config.conf`

At minimum set:

```text
SAMPLESHEET=/home/<user>/Analysis/my_cut_project/config/samplesheet.csv
OUTPUT_DIR=/home/<user>/Analysis/my_cut_project/cutntag2tracks
```

For a human hg38 run, use the validated server references:

```text
INDEX_HG38=/opt/bioinformatics/references/hg38/bowtie2/hg38
FASTA_HG38=/opt/bioinformatics/references/hg38/hg38.fa
CHROM_SIZES_HG38=/opt/bioinformatics/references/hg38/hg38.chrom.sizes
CANONICAL_CONTIGS_HG38=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/hg38/hg38.canonical_contigs.txt
GTF_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf
BLACKLIST_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.blacklist.bed
EFFECTIVE_GENOME_SIZE_HG38=2913022398
TSS_BED_HG38=
CCRE_BED_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz
```

For mouse mm39:

```text
INDEX_MM39=/opt/bioinformatics/references/mm39/bowtie2/mm39
FASTA_MM39=/opt/bioinformatics/references/mm39/mm39.fa
CHROM_SIZES_MM39=/opt/bioinformatics/references/mm39/mm39.chrom.sizes
CANONICAL_CONTIGS_MM39=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/mm39/mm39.canonical_contigs.txt
GTF_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf
BLACKLIST_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.blacklist.bed
EFFECTIVE_GENOME_SIZE_MM39=2654621783
TSS_BED_MM39=
CCRE_BED_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed
```

Remove neither genome block from the template; unused paths are ignored for a
single-genome samplesheet. Review at least peak callers, matched-control policy,
target/control duplicate policy, consensus support, differential modules,
metagene/spike-in choices, cleanup, and resource limits.

For ordinary PE data without experimental spike-in:

```text
SPIKEIN_MODE=none
RUN_METAGENE=false
PEAK_CALLERS=seacr,macs3
TOTAL_CPU_BUDGET=140
RESOURCE_CHECK_MODE=fail
```

For SE data, remove `seacr` from `PEAK_CALLERS`. Add epic2 only when broad or
mixed targets require it and the installed sidecar has been verified.

## 5. Validate metadata without running tools

```bash
CONFIG="$CONFIG_DIR/config.conf"

cutnrun2tracks --config "$CONFIG" --plan
```

Inspect:

```bash
column -t -s $'\t' "$OUTPUT/00_metadata/sample_manifest.tsv" | less -S
column -t -s $'\t' "$OUTPUT/00_metadata/cohort_manifest.tsv" | less -S
column -t -s $'\t' "$OUTPUT/00_metadata/cohort_membership.tsv" | less -S
```

Confirm the reported sequencing-unit, biological-library, and cohort counts and
every target-control assignment.

## 6. Run full preflight

```bash
cutnrun2tracks --config "$CONFIG" --preflight-only
```

Do not continue unless it completes successfully. Review:

```bash
column -t -s $'\t' "$OUTPUT/00_metadata/preflight_status.tsv"
column -t -s $'\t' "$OUTPUT/00_metadata/resource_budget.tsv"
column -t -s $'\t' "$OUTPUT/00_metadata/reference_manifest.tsv" | less -S
column -t -s $'\t' "$OUTPUT/00_metadata/software_versions.tsv" | less -S
```

## 7. Launch the complete workflow

```bash
LOG="$OUTPUT/cutnrun2tracks.nohup.log"
PIDFILE="$OUTPUT/cutnrun2tracks.pid"

nohup cutnrun2tracks --config "$CONFIG" > "$LOG" 2>&1 &
RUN_PID=$!
echo "$RUN_PID" > "$PIDFILE"

echo "PID: $RUN_PID"
echo "Log: $LOG"
```

This single launcher command pins the release and environment internally.

## 8. Monitor progress

```bash
tail -f "$LOG"
```

In another terminal:

```bash
ps -fp "$(cat "$PIDFILE")"
grep '^=== \[' "$OUTPUT/logs/cutnrun2tracks.console.log" | tail -n 20
tail -n 20 "$OUTPUT/00_metadata/workflow_events.tsv"
```

The external `nohup` log and internal console/event logs are intentionally all
retained.

## 9. Verify completion

```bash
grep -E 'COMPLETE|FAILED' "$OUTPUT/00_metadata/workflow_events.tsv" | tail
test -s "$OUTPUT/10_reports/cutnrun2tracks_multiqc_report.html"
test -s "$OUTPUT/10_reports/pipeline_report.html"
test -s "$OUTPUT/00_metadata/final_checksums.sha256"
```

Then inspect the MultiQC report, per-sample peak-call status, consensus summary,
normalization-family status, annotation status, QC plots, and differential
comparison summaries before biological interpretation.

## If the run stops

Do not restart the complete command blindly or delete checkpoints. Follow
[Troubleshooting](16_troubleshooting.md), correct the cause, then resume from
the earliest affected stage:

```bash
RESUME_LOG="$OUTPUT/cutnrun2tracks.resume-<stage>-v0.3.1.log"

nohup cutnrun2tracks \
  --config "$CONFIG" \
  --from-stage <stage> \
  > "$RESUME_LOG" 2>&1 &

echo $! > "$OUTPUT/cutnrun2tracks.resume-<stage>-v0.3.1.pid"
```

`--from-stage` revalidates every earlier checkpoint and reruns the named stage
and all later stages.

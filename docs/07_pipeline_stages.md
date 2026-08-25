# Pipeline stages

`cutnrun2tracks.sh` runs 16 named stages in a fixed order. Each successful
stage receives a signature-and-output JSON checkpoint under
`<OUTPUT_DIR>/.checkpoints/`. The stage names in this page are the accepted
values for `--from-stage` and `--stop-after`.

## Stage contract

| Order | Stage | Main work | Declared output |
|---:|---|---|---|
| 1 | `preflight` | Validate dependencies, configuration, samplesheet, reference compatibility, caller policy, and optional modules | `00_metadata/preflight_status.tsv` |
| 2 | `preprocess` | Merge technical-replicate FASTQs, run raw/trimmed FastQC, optionally trim adapters, and run MultiQC | `02_trimmed_fastq/` |
| 3 | `alignment` | Align to the host reference with Bowtie2; create the optional competitive host-plus-spike alignment | `03_alignment/sorted/` |
| 4 | `filtering` | Mark duplicates once; apply primary/canonical/blacklist and MAPQ filters; expose permissive, intermediate, q30 duplicate-retained, and stringent BAM branches | `03_alignment/analysis/` |
| 5 | `cpm` | Create fragment-aware PE or read-aware SE CPM bedGraph/bigWig tracks and normalization metadata | `04_tracks/cpm/` |
| 6 | `peakcalling` | Run enabled PE SEACR and/or layout-aware MACS3 modes, recording per-sample/caller `SUCCESS`, `EMPTY`, or `ERROR` status under the configured failure policy | `05_peaks/per_sample/` |
| 7 | `consensus` | Build caller/peak-class consensus sets from successful primary peak sets, record exclusions, and enforce biological-sample support | `05_peaks/consensus/` |
| 8 | `spikein` | Count retained spike observations, enforce QC thresholds, and create calibrated host plus spike-control tracks when enabled | `04_tracks/spikein/` |
| 9 | `normalized_tracks` | Build consensus count tables and the enabled DESeq2-consensus and robust-CPM coverage families | `04_tracks/` |
| 10 | `metagene` | Select upstream-normalized tracks and optionally render TSS, TES, and scaled-gene-body profile/heatmap pairs | `06_qc/metagene/status.tsv` |
| 11 | `qc` | Report alignment/complexity, fragment length, FRiP, fingerprints, target-control QC, descriptive TSS signal, and optional experimental ataqv | `06_qc/` |
| 12 | `differential` | Run primary target-only raw-count models and enabled, separately labelled control-aware sensitivity analyses | `08_differential/` |
| 13 | `annotation` | Annotate consensus/differential regions and create UCSC track definitions and an optional IGV session | `07_annotation/`, `09_browser/` |
| 14 | `report` | Assemble the pipeline HTML report and reporting tables | `10_reports/` |
| 15 | `cleanup` | After report success, remove only intermediates selected by the resolved `KEEP_*` policy and record every deletion | `00_metadata/cleanup_status.tsv`, `cleanup_manifest.tsv` |
| 16 | `finalize` | Write final file checksums when enabled and record the terminal checkpoint | `00_metadata/final_checksums.sha256` |

The declared output is the path checked by the stage checkpoint. A stage may
write additional files and directories documented in
[Outputs and recovery](04_outputs_and_recovery.md).

## Optional modules

Optional stages remain in the stage sequence. When a module is disabled, its
script writes an explicit success/skip marker where needed so downstream
checkpoint and reporting behavior stays deterministic.

- `SPIKEIN_MODE=none` disables calibration while retaining the `spikein` stage
  boundary.
- `RUN_METAGENE=false` skips aggregate-signal matrix and plot generation.
- `REQUIRE_ALL_ENABLED_TRACKS=false` records and skips only an unavailable or
  non-normalizable cohort/track family. `true` makes such failures fatal. A
  zero-count sample is never silently removed from its cohort.
- `RUN_DIFFBIND`, `RUN_DESEQ2_ENRICHMENT`,
  `RUN_CONTROL_SUBTRACTED_SENSITIVITY`, and
  `RUN_TARGET_CONTROL_INTERACTION` govern independent differential outputs.
- `RUN_SIMPLE_PEAK_ANNOTATION`, `RUN_CCRE_ANNOTATION`, and
  `WRITE_IGV_SESSION` control optional annotation/browser products.
- `ENABLE_AUTOMATIC_CLEANUP=false` retains intermediates.

## Validation-only runs

Validate the metadata model and write the planned stage order without checking
the existence of FASTQs/references or requiring bioinformatics tools:

```bash
bash cutnrun2tracks.sh \
    --config /absolute/path/to/config.conf \
    --plan
```

Perform the complete dependency/reference/input audit and stop before
preprocessing:

```bash
bash cutnrun2tracks.sh \
    --config /absolute/path/to/config.conf \
    --preflight-only
```

## Controlled reruns

Reuse hash-validated outputs before one stage, then force that stage and all
later stages:

```bash
bash cutnrun2tracks.sh \
    --config /absolute/path/to/config.conf \
    --from-stage normalized_tracks
```

Earlier checkpoint signatures are adopted only after their recorded output
sizes and SHA-256 hashes pass validation. Adoption is recorded in each JSON
checkpoint. Missing or changed earlier outputs stop the run. Use this explicit
override only when the change cannot affect stages before the selected boundary.

Stop after a named stage:

```bash
bash cutnrun2tracks.sh \
    --config /absolute/path/to/config.conf \
    --stop-after qc
```

Do not manually copy or edit a checkpoint between runs. Its signature covers the
workflow version, scripts, shared modules, resolved configuration, sanitized
samplesheet, and reference manifest, and its output hashes are specific to the
run directory.

## Scientific boundaries between stages

- `cpm`, `spikein`, and `normalized_tracks` produce visualization tracks.
  `differential` uses raw fragment/read counts, not bigWig values.
- `normalized_tracks` writes `normalized_track_family_status.tsv` and
  per-sample consensus-count diagnostics. An upstream cohort skip is propagated
  as a differential skip, allowing later QC and reporting to run.
- Matched controls are resolved before processing and are used for background
  peak calling and QC. The primary differential analysis is target-only.
- Technical replicates merge during `preprocess`; biological replicates remain
  distinct through consensus support and differential modeling.
- `metagene` consumes completed bigWigs and never renormalizes them.
- `annotation` genome-sorts both consensus and GTF-derived BED records using
  the configured chromosome-sizes order before nearest-gene lookup.
- `cleanup` cannot start unless the final report exists, and `finalize` records
  the post-cleanup deliverable tree.

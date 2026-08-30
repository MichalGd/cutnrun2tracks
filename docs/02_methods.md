# Methods overview

[Documentation index](README.md) | [Pipeline stages](07_pipeline_stages.md)

This page is a map of the analytical method implemented by
`cutnrun2tracks` 0.3.1. Detailed, implementation-specific descriptions are in
the linked pages. The resolved configuration and metadata under
`00_metadata/` remain authoritative for an individual run.

## Analysis unit

Rows with the same `sample_id` and biological `replicate` are technical
sequencing units of one biological library. Their gzip FASTQs are concatenated
before trimming; they are not treated as independent replicates. Biological
replicates remain separate throughout alignment, peak calling, consensus
construction, QC, and differential analysis.

For paired-end (PE) libraries, the intended signal unit is one retained proper
fragment, represented by the outer aligned coordinates and counted once for
coverage and interval QC. Single-end (SE) signal is one retained primary read;
the workflow does not infer nucleosomal fragments from SE data.

## Processing outline

1. Validate configuration, samples, controls, references, tools, and the
   declared CPU budget.
2. Run FastQC on technical FASTQs, merge technical units, optionally trim with
   Trim Galore/cutadapt, and run FastQC on merged raw and trimmed libraries.
3. Align to the host or competitive host-plus-spike reference with Bowtie2.
4. Mark duplicates and construct four canonical-contig, blacklist-filtered BAM
   branches that vary MAPQ and duplicate retention.
5. Select the q30 branch requested by each sample's duplicate policy as its
   `analysis` BAM.
6. Generate CPM tracks, per-sample peak calls, biological-support consensus
   regions, cohort-normalized tracks, and optional spike-in tracks.
7. Generate sequencing and assay-performance QC, metagene plots, differential
   analyses, genomic annotation, browser files, and final reports.

## Filtering branches

| Name used in the documentation | MAPQ | Duplicates | Main purpose |
|---|---:|---|---|
| permissive | 0 | retained | sensitivity to ambiguous mappings and duplicates |
| intermediate | 0 | removed | sensitivity to ambiguous mappings after deduplication |
| q30 duplicate-retained | 30 | retained | analysis option for sparse/low-input target libraries |
| stringent | 30 | removed | conservative analysis and normalization option |

The samplesheet `duplicate_policy` selects the q30 retained or q30 removed BAM
as the analysis BAM. Target and control defaults are independently configurable.
All four branches are made from the marked BAM rather than serially from one
another. See [Reference filtering and blacklists](10_references_blacklist_and_filtering.md).

## Detailed method pages

- [Peak calling and consensus](09_peak_calling.md): SEACR, MACS3, epic2,
  matched controls, defaults, caller failure isolation, and biological support.
- [Reference filtering and blacklists](10_references_blacklist_and_filtering.md):
  exact server files, canonical-contig selection, MAPQ/flag/duplicate filters,
  and blacklist mechanics.
- [Quality control](11_quality_control.md): separate sequencing-QC and
  assay-performance-QC inventories covering every implemented metric.
- [Tracks and spike-in normalization](12_tracks_and_normalization.md): every
  formal track family, formulas, intended comparisons, and calibration QC.
- [Differential binding](03_differential_enrichment.md): target-only DESeq2,
  DiffBind, blocking, control-subtracted and interaction sensitivities.
- [Genomic annotation](13_genomic_annotation.md): consensus lookup,
  differential-table propagation, and mutually exclusive feature summaries for
  every successful peak set.
- [Metagene plots](06_metagene.md): TSS, TES, and scaled gene-body aggregate
  profiles.

## Interpretation principle

The workflow deliberately reports multiple callers, filtering policies,
normalizations, and QC measurements. They are complementary sensitivity views,
not interchangeable duplicates of one result. No single QC number or caller is
promoted to a universal CUT&RUN/CUT&Tag pass/fail rule. The analyst should
record a project-specific inclusion decision using the combined sequencing,
assay, replicate, control, and biological evidence.

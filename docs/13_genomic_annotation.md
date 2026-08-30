# Genomic annotation

[Documentation index](README.md) | [Differential binding](03_differential_enrichment.md)

Annotation runs after peak calling and statistical testing. It adds descriptive
genomic context but never changes peak coordinates, consensus support, raw
counts, fold changes, P values, or adjusted P values. Version 0.3.1 implements
two complementary layers:

1. lightweight nearest-gene and cCRE lookup for the primary consensus and
   gzipped differential tables; and
2. exhaustive mutually exclusive feature assignment and composition summaries
   for every successful per-sample peak set from every enabled caller, plus the
   primary consensus when requested.

## Reference files

Annotation resolves references by each samplesheet row's genome:

```text
GTF_HG38, GTF_MM39
CHROM_SIZES_HG38, CHROM_SIZES_MM39
CCRE_BED_HG38, CCRE_BED_MM39
```

The documented shared-server deployment uses:

| Genome | GTF | cCRE |
|---|---|---|
| hg38 | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf` | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz` |
| mm39 | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed` |

The exact run-specific paths and checksums are in
`00_metadata/reference_manifest.tsv`. Record the GTF and cCRE source releases;
filenames alone do not identify annotation provenance. Human/mouse GTF concepts
are described by the [GENCODE project](https://www.gencodegenes.org/), and the
ENCODE registry of candidate cis-regulatory elements is described by the ENCODE
Project Consortium ([Nature 2020, 10.1038/s41586-020-2493-4](https://doi.org/10.1038/s41586-020-2493-4)).

## Layer 1: primary-consensus lookup

For each cohort with a primary consensus, `RUN_SIMPLE_PEAK_ANNOTATION=true`
does the following:

1. extract GTF rows whose feature type is exactly `gene`;
2. convert their 1-based GTF spans to zero-based BED, retaining gene name,
   gene ID, and strand;
3. sort genes and consensus peaks in configured chromosome-size order; and
4. run `bedtools closest -d`.

The result is:

```text
07_annotation/<cohort>/consensus/<cohort>.nearest_gene.tsv
```

Distance is unsigned and zero when a consensus interval overlaps the whole gene
span. This is nearest **gene span**, not nearest TSS, and proximity does not
establish regulatory targeting.

With `RUN_CCRE_ANNOTATION=true`, `bedtools intersect -wao` preserves every
primary-consensus/cCRE overlap, all supplied cCRE BED fields, and overlap length:

```text
07_annotation/<cohort>/consensus/<cohort>.ccre_reference_overlaps.tsv
```

Non-overlapping consensus intervals remain present with BEDTools placeholder
fields. BEDTools interval algorithms are described by Quinlan and Hall
([Bioinformatics 2010, 10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)).

### Differential-table propagation

After consensus lookup, the workflow scans gzipped TSV files below the cohort's
`08_differential/` tree. Any table containing `region_id` receives an adjacent
`*.annotated.tsv.gz` copy with:

| Added column | Definition |
|---|---|
| `nearest_gene_name` | name of nearest whole-gene interval |
| `nearest_gene_id` | its GTF gene identifier |
| `distance_to_gene` | unsigned distance to gene span; zero for overlap |
| `ccre_reference_overlaps` | unique overlapping values from cCRE BED column 4 |

Original result tables remain the statistical source of truth. Direct DESeq2
and interaction results are gzipped and are annotated. Current DiffBind contrast
TSVs are uncompressed and therefore do not receive automatic copies in v0.3.1.

## Layer 2: all-peak genomic-feature composition

With `RUN_FEATURE_ANNOTATION_SUMMARY=true`, the workflow reads
`05_peaks/per_sample/<sample>/caller_status.tsv` and annotates every valid peak
file marked `SUCCESS`, including non-primary caller/class sensitivity results.
Controls are not peak-called and are absent. With
`PEAK_ANNOTATION_INCLUDE_CONSENSUS=true`, successful primary consensus BEDs are
added as separate `entity_type=consensus` entries.

### Feature construction

GTF `gene` and `exon` rows are converted to these strand-aware features:

- **promoter**: by default 2,000 bp upstream and 500 bp downstream of TSS,
  controlled by `PEAK_ANNOTATION_PROMOTER_UPSTREAM` and `_DOWNSTREAM`;
- **exon**: merged exon intervals for each gene;
- **intron**: portions of the GTF gene span not covered by its merged exons;
- **gene_end**: the 2,000 bp interval immediately downstream of TES, controlled
  by `PEAK_ANNOTATION_GENE_END_WINDOW`.

cCRE rows are classified heuristically from their supplied text fields:

- labels containing `PLS` or `promoter` become **promoter**;
- labels containing `dELS`, `pELS`, `enhancer`, or ` ELS` become **enhancer**;
- remaining cCRE records become **other_regulatory**.

This mapping depends on the cCRE BED vocabulary and is not a de novo enhancer
prediction. A cCRE overlap also does not prove activity in the assayed cell type
or establish an enhancer-gene link.

### Exclusive primary assignment

A peak can overlap several genes/features. The all-overlaps table preserves all
such relationships, while the summary assigns exactly one primary category.
The default precedence is:

```text
promoter > enhancer > exon > intron > gene_end >
other_regulatory > intergenic > unclassified
```

`PEAK_ANNOTATION_FEATURE_PRECEDENCE` may reorder the same complete category
set. The algorithm first selects the highest-precedence category, then the
greatest overlap in base pairs, then a deterministic feature ID tie-break.
Peaks with no feature overlap are `intergenic`; peaks on contigs absent from the
configured chromosome sizes are `unclassified`.

Every peak also receives the nearest GTF gene to its midpoint and a signed,
transcription-oriented TSS distance: negative is upstream and positive is
downstream. This nearest-TSS measure is distinct from Layer 1's unsigned nearest
whole-gene distance.

### Tables and plots

Outputs under `07_annotation/feature_summary/` are:

| File | Contents |
|---|---|
| `peak_feature_assignments.tsv.gz` | one row per peak with primary category, all overlapping categories, primary feature/gene, overlap bp, nearest gene, and signed TSS distance |
| `peak_feature_all_overlaps.tsv.gz` | one row per peak-feature overlap; preserves nonexclusive evidence |
| `peak_feature_summary.tsv` | combined count, fraction, percentage, assigned peak bp, bp fraction, totals, and color |
| `peak_feature_counts.tsv` | absolute number of peaks in every category |
| `peak_feature_fractions.tsv` | fraction and percentage of peaks in every category |
| `peak_feature_bp_coverage.tsv` | total peak width and fraction assigned to every category |
| `peak_annotation_status.tsv` | source peak file, valid/invalid/unclassified counts, enhancer-reference status, and completion state |
| `peak_feature_colors.tsv` | stable category/color mapping |

For every caller/peak-class combination, the workflow creates
`peak_feature_composition.<caller>.<class>.{png,pdf,svg}` by default. Each sample
is one horizontal color-coded stacked bar in three panels: absolute peak count,
fraction of peaks, and fraction of peak-covered base pairs. The bp summary
assigns each peak's entire width to its exclusive primary category; it is not a
sum of the potentially overlapping bases in the all-overlaps table.

## Configuration

```text
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=true
RUN_FEATURE_ANNOTATION_SUMMARY=true
PEAK_ANNOTATION_PROMOTER_UPSTREAM=2000
PEAK_ANNOTATION_PROMOTER_DOWNSTREAM=500
PEAK_ANNOTATION_GENE_END_WINDOW=2000
PEAK_ANNOTATION_FEATURE_PRECEDENCE=promoter,enhancer,exon,intron,gene_end,other_regulatory,intergenic,unclassified
PEAK_ANNOTATION_PLOT_FORMATS=png,pdf,svg
PEAK_ANNOTATION_INCLUDE_CONSENSUS=true
RUN_MOTIF_ENRICHMENT=false
```

Motif enrichment is not implemented; enabling it is rejected. If cCRE is
disabled or unavailable, enhancer annotation is recorded as `not_evaluated`
rather than silently calling every non-genic peak intergenic.

## Interpretation limitations

- “Nearest” is a geometric relationship, not a validated target gene.
- GTF completeness, transcript/gene definitions, and cCRE release directly
  affect category counts.
- Exclusive precedence is necessary for a stacked composition plot but hides
  secondary overlaps; consult `peak_feature_all_overlaps.tsv.gz` when biology
  depends on overlapping labels.
- The feature summary describes detected intervals, not their statistical
  confidence, signal amplitude, or differential status.
- Per-sample callers with different peak widths or thresholds can have different
  category compositions for technical as well as biological reasons.
- No chromatin-state segmentation, enhancer-to-gene linking, ontology analysis,
  or motif analysis is inferred by this stage.

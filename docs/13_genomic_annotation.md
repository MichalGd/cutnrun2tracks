# Genomic annotation

[Documentation index](README.md) | [All reference families](17_reference_preparation_and_provenance.md) | [Differential binding](03_differential_enrichment.md)

`cutnrun2tracks` adds genomic context after peak calling and statistical
testing. Annotation never changes peak coordinates, consensus support, raw
counts, fold changes, P values, or adjusted P values. It has two complementary
layers:

1. a lightweight nearest-gene and cCRE lookup for each primary consensus, with
   propagation into compatible differential tables; and
2. exhaustive, mutually exclusive feature assignment and composition summaries
   for every successful per-sample peak set from every enabled caller, plus the
   primary consensus when requested.

The second layer is broader than the corresponding ATACseq2tracks summary: it
annotates per-sample MACS3, SEACR, and epic2 results, not only a differential
region universe. Its column names and tie-breaking rules are therefore
documented here rather than assumed to be identical to ATACseq2tracks.

## Scope and propagation

With the default annotation options, the workflow considers:

- every valid target peak file whose `caller_status.tsv` row is `SUCCESS`,
  including non-primary caller and peak-class sensitivity results;
- every successful primary consensus when
  `PEAK_ANNOTATION_INCLUDE_CONSENSUS=true`;
- every interval in those files, not only peaks passing a later differential
  threshold; and
- primary-consensus annotations joined by stable `region_id` to gzipped
  differential TSVs that contain that identifier.

Controls are not peak-called and therefore do not appear as per-sample peak
entities. Different callers and peak classes define different interval
universes; rows should not be assumed to correspond one-to-one across MACS3
narrow, MACS3 broad, SEACR, epic2, or consensus results.

## Genome-specific sources and construction

### Deployed paths and upstream identities

The shared-server configuration resolves annotation references by the genome in
each samplesheet row:

| Genome | Config key | Deployed file used by the workflow | Intended upstream identity |
|---|---|---|---|
| hg38 | `GTF_HG38` | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf` | GENCODE release 42 comprehensive primary-assembly GTF, `gencode.v42.primary_assembly.annotation.gtf.gz`, GRCh38.p13 |
| hg38 | `CCRE_BED_HG38` | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz` | ENCODE4 expanded Registry, `Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz`, native GRCh38 |
| mm39 | `GTF_MM39` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf` | GENCODE release M31 comprehensive primary-assembly GTF, `gencode.vM31.primary_assembly.annotation.gtf.gz`, GRCm39 |
| mm39 | `CCRE_BED_MM39` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed` | UCSC/ENCODE3 mm10 `encodeCcreCombined.bb`, converted to BED and lifted to mm39 |

The GENCODE files are available from the official [human release 42](https://www.gencodegenes.org/human/release_42.html)
and [mouse release M31](https://www.gencodegenes.org/mouse/release_M31.html)
pages and the corresponding EMBL-EBI FTP directories. The installed
`annotation.gtf` files are decompressed aliases of the named primary-assembly
GTFs; no liftOver or coordinate conversion is applied to either GTF. GTF
coordinates are converted from one-based closed intervals to zero-based,
half-open intervals only in memory during annotation.

The alias name alone is not proof of source identity. The shared deployment
must retain its construction manifest, and every run records the resolved path,
byte size, and SHA-256 digest in:

```text
00_metadata/reference_manifest.tsv
```

That run manifest proves which bytes were used. It does not reconstruct the
download URL or release label, so the server-level provenance record must be
retained with the reference files.

### Human hg38 construction

#### GTF

1. Download
   [`gencode.v42.primary_assembly.annotation.gtf.gz`](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz)
   from GENCODE release 42.
2. Verify the downloaded MD5 against the release directory's `MD5SUMS` file.
3. Decompress it without modifying records and install the result as
   `references/hg38/annotation.gtf`.
4. Record both the upstream MD5 and installed SHA-256, byte size, source URL,
   release, assembly, and installation date in the shared reference manifest.
5. Confirm that `gene` and `exon` records have `gene_id`, and that contig names
   agree with `CHROM_SIZES_HG38`.

GENCODE reports 62,696 genes for release 42 on the main chromosomes. This is a
release benchmark, not an exact raw-line count for the primary-assembly GTF,
which can also contain scaffold records. Acceptance of the installed file must
therefore use counts computed from that exact file rather than substituting the
summary statistic.

#### cCRE

The human cCRE source is the native GRCh38 ENCODE4 expanded-registry file
[`Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz`](https://users.moore-lab.org/ENCODE-cCREs/Supplementary-Data/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz).
No liftOver is performed. The ATACseq2tracks preparation utility
`utilities/prepare_encode4_hg38_ccre.sh`, used to construct the shared
reference, performs the following acceptance checks:

- a valid gzip stream;
- exactly **2,348,854** non-comment records in this canonical download;
- canonical GRCh38 chromosome names and valid zero-based BED coordinates;
- non-empty accessions and recognized ENCODE4 class labels;
- unchanged SHA-256 before and after installation; and
- a provenance TSV containing source URL, filename, download time, assembly,
  registry name, record count, BED columns, and SHA-256.

The deployed `hg38.ccre.bed.gz` should be a byte-preserving alias of that
validated download. A broader registry total from a paper or website must not
replace the file-level count above.

### Mouse mm39 construction

#### GTF

1. Download
   [`gencode.vM31.primary_assembly.annotation.gtf.gz`](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz)
   from GENCODE release M31.
2. Verify it against the release directory's `MD5SUMS` file.
3. Decompress it without changing records and install it as
   `references/mm39/annotation.gtf`.
4. Record upstream and installed checksums, byte size, URL, release, assembly,
   and installation date; validate `gene_id` and chromosome compatibility.

GENCODE reports 56,923 genes for release M31 on the main chromosomes. As for
human, this is a release benchmark rather than the exact primary-assembly GTF
line count.

#### cCRE

The deployed mouse cCRE file is not a native mm39 ENCODE4 registry. It was
constructed from the archival ENCODE3/UCSC mm10 combined track:

1. download `encodeCcreCombined.bb` from
   `http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/`;
2. convert bigBed to BED with UCSC Kent `bigBedToBed`;
3. download the UCSC `mm10ToMm39.over.chain.gz` chain;
4. map the BED with `liftOver -bedPlus=9 -tab`, preserving the first nine BED
   fields plus all custom fields;
5. retain rejected intervals in a separate unmapped file; and
6. coordinate-sort mapped intervals with `LC_COLLATE=C sort -k1,1 -k2,2n` to
   create the installed `mm39.ccre.bed`.

The construction is recorded in the historical utility
[`creating_ENCODE_cCRES_mm39_bigBed_track.sh`](https://github.com/MichalGd/ATAC-seq/blob/main/utilitiies/creating_ENCODE_cCRES_mm39_bigBed_track.sh).
The mm10 source contains **339,815** cCREs and uses the ENCODE3 classes `PLS`,
`pELS`, `dELS`, `DNase-H3K4me3`, and `CTCF-only`. The exact installed mm39
record count is the number that mapped successfully, not 339,815; mapped and
unmapped source records should reconcile to the source total. The installed
count, unmapped count, chain checksum, Kent-tool versions, BED checksum, and
construction date must remain in the server provenance record.

LiftOver can reject, split, or reposition intervals. Human and mouse therefore
use different registry generations and construction routes, and their cCRE
class counts are not directly comparable.

### Reference acceptance audit

On the server, the following read-only audit captures the exact deployed
counts and digests. Store its output beside the shared reference manifest:

```bash
HG_GTF=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf
MM_GTF=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf
HG_CCRE=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz
MM_CCRE=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed

sha256sum "$HG_GTF" "$MM_GTF" "$HG_CCRE" "$MM_CCRE"

for gtf in "$HG_GTF" "$MM_GTF"; do
    printf '%s\tgenes=' "$gtf"
    awk -F '\t' '$0 !~ /^#/ && $3 == "gene" {n++} END {print n+0}' "$gtf"
    printf '%s\texons=' "$gtf"
    awk -F '\t' '$0 !~ /^#/ && $3 == "exon" {n++} END {print n+0}' "$gtf"
done

printf '%s\trecords=' "$HG_CCRE"
gzip -cd "$HG_CCRE" | awk '$0 !~ /^(#|track|browser)/ && NF {n++} END {print n+0}'
printf '%s\trecords=' "$MM_CCRE"
awk '$0 !~ /^(#|track|browser)/ && NF {n++} END {print n+0}' "$MM_CCRE"
```

Expected human cCRE output is 2,348,854. GTF and lifted-mouse counts are
deployment-specific acceptance values and should be compared with the retained
construction manifest. Preflight rejects missing or empty references, but it
does not independently certify that an arbitrary configured BED is the claimed
registry release.

## GTF-derived feature construction

The exhaustive layer imports only GTF rows whose feature type is exactly
`gene` or `exon` and that contain `gene_id`. `gene_name` is used when present,
otherwise `gene_id` becomes the display name. Peak strands are ignored because
CUT&RUN/CUT&Tag peak intervals are treated as unstranded.

For every gene, the workflow constructs:

| Feature category | Implemented definition |
|---|---|
| `promoter` | Strand-aware window from 2,000 bp upstream through 500 bp downstream of the GTF gene TSS by default |
| `exon` | Union of all GTF exon intervals belonging to the gene; overlapping exons are merged |
| `intron` | Parts of the GTF gene span not covered by the merged exon union |
| `gene_end` | Strand-aware interval immediately downstream of the GTF gene TES, 2,000 bp by default |

The default promoter definition is controlled by:

```bash
PEAK_ANNOTATION_PROMOTER_UPSTREAM=2000
PEAK_ANNOTATION_PROMOTER_DOWNSTREAM=500
PEAK_ANNOTATION_GENE_END_WINDOW=2000
```

For a plus-strand gene, the promoter is
`[TSS - upstream, TSS + downstream)` and the gene-end interval starts at the
TES. For a minus-strand gene, the orientation is reversed: the promoter extends
`downstream` bases toward lower genomic coordinates and `upstream` bases toward
higher coordinates, while the gene-end interval lies below the gene start.
Coordinates are clipped at zero.

Because features are constructed per gene, a peak can overlap promoters,
exons, introns, or gene ends from different overlapping genes simultaneously.
The all-overlaps output preserves these relationships; the primary category is
only a deterministic summary.

## cCRE regulatory classification

The cCRE registry has already classified candidate regulatory elements from
integrated biochemical evidence. `cutnrun2tracks` does not infer chromatin
states from the current CUT&RUN/CUT&Tag samples. It searches the text in BED
columns 4 onward case-insensitively and maps source labels to its broader plot
categories:

| Source label | Registry interpretation | `cutnrun2tracks` category |
|---|---|---|
| `PLS` | promoter-like signature | `promoter` |
| `pELS` | proximal enhancer-like signature | `enhancer` |
| `dELS` | distal enhancer-like signature | `enhancer` |
| `CA-H3K4me3` | accessible, H3K4me3-associated ENCODE4 class | `other_regulatory` |
| `CA-CTCF` | accessible and CTCF-associated ENCODE4 class | `other_regulatory` |
| `CA-TF` | accessible, TF-bound ENCODE4 class | `other_regulatory` |
| `CA` | accessible ENCODE4 class lacking the defining marks above | `other_regulatory` |
| `TF` | TF-bound ENCODE4 class | `other_regulatory` |
| `DNase-H3K4me3` | legacy ENCODE3 H3K4me3-associated non-PLS class | `other_regulatory` |
| `CTCF-only` | legacy ENCODE3 CTCF-associated class | `other_regulatory` |
| `CTCF-bound` | additional CTCF-bound tag when present | `other_regulatory` |
| unrecognized cCRE text | supplied regulatory record without a recognized promoter/enhancer token | `other_regulatory` |

The parser also recognizes literal `promoter`, `enhancer`, and ELS-like text.
The feature identifier is BED column 4, or a generated `ccre:<row-number>` ID
when that field is empty.

Unlike ATACseq2tracks, the exhaustive CUT annotation does not emit separate
`ccre_primary_class`, `ccre_all_classes`, or `enhancer_like` columns. It maps
cCREs into the categories above and preserves each cCRE as a row in
`peak_feature_all_overlaps.tsv.gz`. The lightweight consensus layer separately
preserves the complete supplied BED record in
`*.ccre_reference_overlaps.tsv`.

## Exclusive primary classification and tie-breaking

Each valid peak receives exactly one `primary_category`. The default precedence
is:

```text
promoter > enhancer > exon > intron > gene_end >
other_regulatory > intergenic > unclassified
```

Selection is deterministic:

1. collect every GTF-derived or cCRE feature with at least one overlapping
   base;
2. select the highest-precedence category present;
3. within that category, select the feature with the greatest overlap in base
   pairs; and
4. if overlap is still tied, select the lexicographically smallest feature ID.

This differs deliberately from ATACseq2tracks' primary-cCRE selection, which
first maximizes overlap and then applies cCRE-class priority. In
`cutnrun2tracks`, category precedence is applied before overlap because GTF and
cCRE features share one exclusive composition system. A cCRE `PLS` and a GTF
promoter both enter the `promoter` category and are then resolved by overlap and
feature ID.

Peaks with no feature overlap are `intergenic`. Peaks on contigs absent from
the configured chromosome-sizes file are `unclassified`, even if a supplied
reference happens to contain a feature on that contig. The precedence can be
reordered with `PEAK_ANNOTATION_FEATURE_PRECEDENCE`, but the validator requires
the same complete category set with no duplicates.

## Nearest-gene definitions

The two annotation layers answer different questions:

| Layer | Reference point | Distance |
|---|---|---|
| Lightweight primary-consensus lookup | nearest whole GTF `gene` interval via `bedtools closest -d` | unsigned; zero for gene-span overlap |
| Exhaustive feature assignment | peak midpoint to nearest strand-aware GTF gene TSS | signed in transcriptional orientation; negative upstream, positive downstream |

The exhaustive nearest-TSS tie is resolved deterministically by absolute
distance, then gene ID and gene name. Neither definition establishes a
regulatory target.

## Output columns

### Exhaustive exclusive assignments

`07_annotation/feature_summary/peak_feature_assignments.tsv.gz` contains one
row per valid peak:

| Column | Definition |
|---|---|
| `entity_type` | `sample` or `consensus` |
| `entity_id` | sample key or cohort identifier |
| `sample_key` | biological sample key; `.` for consensus entities |
| `cohort_id` | resolved target cohort |
| `genome`, `factor`, `condition`, `replicate` | samplesheet-derived identity fields |
| `caller`, `peak_class` | originating caller and narrow/broad/stringent class |
| `peak_file` | exact annotated source file |
| `peak_id` | source BED column 4 or generated stable row ID |
| `chrom`, `start`, `end`, `width` | zero-based half-open peak coordinates and width |
| `primary_category` | single precedence-based category |
| `overlapping_categories` | all unique overlapping categories, ordered by configured precedence |
| `primary_feature_id` | winning GTF-derived feature ID or cCRE BED-column-4 ID |
| `primary_gene_id`, `primary_gene_name` | gene attached to the winning GTF feature; `.` for cCRE/intergenic assignments |
| `primary_overlap_bp` | overlap length with the winning feature |
| `nearest_gene_id`, `nearest_gene_name` | nearest gene by midpoint-to-TSS distance |
| `nearest_tss_signed_distance` | signed transcription-oriented midpoint-to-TSS distance in bp |
| `enhancer_annotation` | whether a cCRE reference was available (`evaluated` or the recorded unavailable state) |

### Nonexclusive all-overlaps table

`peak_feature_all_overlaps.tsv.gz` has one row for every peak-feature overlap.
It repeats entity and peak coordinates and adds:

| Column | Definition |
|---|---|
| `feature_category` | mapped promoter/enhancer/exon/intron/gene-end/other-regulatory category |
| `feature_id` | GTF-derived feature ID or cCRE identifier |
| `feature_start`, `feature_end` | feature coordinates |
| `overlap_bp` | intersection length |
| `gene_id`, `gene_name`, `strand` | GTF gene metadata, or `.` for cCRE features |

This table is the authoritative source when secondary overlaps matter.

### Composition tables and plots

| File | Contents |
|---|---|
| `peak_feature_summary.tsv` | counts, fractions, percentages, assigned peak bp, bp fractions, totals, and colors |
| `peak_feature_counts.tsv` | absolute peak counts for every category and entity |
| `peak_feature_fractions.tsv` | count fractions and percentages |
| `peak_feature_bp_coverage.tsv` | total peak width assigned to each exclusive category and its fraction |
| `peak_annotation_status.tsv` | input file, valid/invalid/unclassified counts, cCRE evaluation state, status, and reason |
| `peak_feature_colors.tsv` | stable category-to-color mapping |

For each caller/peak-class combination, the workflow writes
`peak_feature_composition.<caller>.<class>.{png,pdf,svg}`. Each sample is one
horizontal color-coded stacked bar in three panels: absolute peak count,
fraction of peaks, and fraction of peak-covered base pairs. The bp panel assigns
each peak's entire width to its exclusive primary category; it does not sum
overlap lengths from the nonexclusive table.

### Lightweight consensus and differential outputs

For each primary consensus, Layer 1 creates:

```text
07_annotation/<cohort>/consensus/<cohort>.nearest_gene.tsv
07_annotation/<cohort>/consensus/<cohort>.ccre_reference_overlaps.tsv
```

The cCRE file is produced with `bedtools intersect -wao`, preserves every
consensus interval, all supplied cCRE BED fields, and overlap length. The
nearest-gene file is produced after both genes and consensus intervals are
sorted in configured chromosome-size order.

Compatible gzipped differential TSVs containing `region_id` receive adjacent
`*.annotated.tsv.gz` copies with:

| Added column | Definition |
|---|---|
| `nearest_gene_name` | name of nearest whole-gene interval |
| `nearest_gene_id` | GTF identifier of that interval |
| `distance_to_gene` | unsigned distance to gene span; zero for overlap |
| `ccre_reference_overlaps` | unique values from cCRE BED column 4 |

Original differential tables remain the statistical source of truth. In
version 0.3.1, uncompressed DiffBind contrast TSVs do not receive automatic
annotated copies.

## Configuration and opt-out behavior

```bash
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

With cCRE annotation enabled, preflight requires a non-empty genome-matched
reference. Generic preflight does not reproduce the stricter release-specific
construction audit above. If cCRE annotation is disabled or unavailable,
GTF-derived promoter/exon/intron/gene-end classification can still run and the
enhancer evaluation state is recorded rather than silently treating all
non-genic peaks as biologically intergenic. Motif enrichment is not
implemented; enabling it is rejected.

## Interpretation limits

- A coordinate overlap is descriptive evidence, not proof that a reference
  element is active in the assayed cell type.
- `enhancer` means overlap with a supplied pELS/dELS or enhancer-labelled
  reference interval. It is not a functional enhancer assay and does not link
  the element to a target gene.
- A nearest gene is a geometric neighbour, not a validated regulatory target.
  Long-range regulation, chromatin contacts, and intervening genes are not
  modelled.
- Cell-type-agnostic cCRE registries integrate evidence across many biosamples;
  an element can be inactive in the current experiment.
- The lifted mouse registry and native human registry use different releases,
  class vocabularies, and construction histories. Cross-species category
  fractions are not directly equivalent.
- Broad peaks can overlap many genes and regulatory elements. Consult the
  all-overlaps table and original intervals rather than only the primary label.
- Exclusive precedence is necessary for stacked summaries but hides secondary
  relationships. Reordering it changes composition counts without changing the
  underlying peak calls.
- Category composition depends on GTF release, cCRE release, assembly,
  chromosome set, peak caller, peak width, and thresholds. It does not by itself
  measure statistical confidence, signal intensity, or differential binding.
- No chromatin-state segmentation, enhancer-to-gene linking, ontology analysis,
  motif analysis, or causal regulatory inference is performed.

## References

- Frankish A, et al. [GENCODE: reference annotation for the human and mouse
  genomes in 2023](https://doi.org/10.1093/nar/gkac1071). *Nucleic Acids
  Research* (2023).
- ENCODE Project Consortium. [Expanded encyclopaedias of DNA elements in the
  human and mouse genomes](https://doi.org/10.1038/s41586-020-2493-4).
  *Nature* (2020).
- Moore JE, Pratt HE, Fan K, et al. [An expanded registry of candidate
  cis-regulatory elements](https://doi.org/10.1038/s41586-025-09909-9).
  *Nature* (2026).
- ENCODE Project. [Candidate Cis-Regulatory Elements pipeline](https://www.encodeproject.org/pipelines/ENCPL751FOQ/)
  and [cCRE subtype glossary](https://www.encodeproject.org/glossary/).
- Kent WJ, et al. [The Human Genome Browser at UCSC](https://doi.org/10.1101/gr.229102).
  *Genome Research* (2002). See also the UCSC Kent utilities used for bigBed
  conversion and liftOver.
- Quinlan AR, Hall IM. [BEDTools: a flexible suite of utilities for comparing
  genomic features](https://doi.org/10.1093/bioinformatics/btq033).
  *Bioinformatics* (2010).

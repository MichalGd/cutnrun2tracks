# Metagene aggregate-signal module

The metagene stage creates a profile and matching heatmap for every selected
sample, gene set, and mode. Supported modes are TSS-centered, TES-centered, and
scaled gene body. It uses deepTools `computeMatrix`, `plotProfile`, and
`plotHeatmap` with BED12 exon blocks and `--metagene`.

The stage consumes final bigWigs. It does not calculate CPM, RPGC, spike-in, or
other normalization. The chosen upstream normalization is recorded in every
task name and metadata row.

## Enable the workflow stage

Prepare references first, create a `gene_sets.tsv` manifest, and set:

```bash
RUN_METAGENE=true
METAGENE_GENE_SET_MANIFEST=/refs/metagene/gene_sets.tsv
METAGENE_GENE_SETS=protein_coding,broadly_expressed
METAGENE_MODES=tss,tes,gene_body
METAGENE_TRACK_FAMILY=auto
```

With `METAGENE_TRACK_FAMILY=auto`, spike-in-normalized host tracks are selected
when spike-in is enabled; otherwise CPM tracks are selected. A missing or failed
spike-in track is fatal unless `METAGENE_ALLOW_CPM_FALLBACK=true` is explicitly
set. Such fallbacks retain the `CPM` normalization label.

Controls are excluded by default. Set `METAGENE_INCLUDE_CONTROLS=true` to plot
them.

## Protein-coding BED12 reference

```bash
python common/metagene/prepare_metagene_reference.py \
  --gtf /refs/gencode.annotation.gtf.gz \
  --assembly GRCh38 \
  --annotation-source gencode \
  --annotation-release vNN \
  --chrom-sizes /refs/GRCh38.chrom.sizes \
  --canonical-contigs /refs/GRCh38.canonical_contigs.txt \
  --blacklist /refs/GRCh38.blacklist.bed \
  --transcript-policy canonical_then_longest \
  --min-gene-span 1000 \
  --min-spliced-length 500 \
  --blacklist-policy gene_span \
  --overlap-policy keep \
  --output-dir /refs/metagene/GRCh38/gencode_vNN
```

The builder selects one protein-coding transcript per gene, using MANE Select,
Ensembl Canonical, APPRIS, CDS length, spliced length, and stable transcript ID
as deterministic priorities. It emits unfiltered and filtered BED12 files,
gene metadata, exclusion reasons, and a checksummed reference manifest.

To remove overlapping or nearby loci, use `--overlap-policy exclude` and
`--adjacency-bp N`. Collisions are evaluated against the complete eligible
protein-coding universe.

## HPA reference subset

```bash
python common/metagene/build_hpa_reference_subset.py \
  --hpa-table /refs/HPA/proteinatlas.tsv.zip \
  --hpa-release 25.1 \
  --human-gene-model /refs/metagene/GRCh38/gencode_vNN/protein_coding.filtered.bed12 \
  --gene-set-id broadly_expressed \
  --specificity "Low tissue specificity" \
  --distribution "Detected in all" \
  --output-dir /refs/metagene/HPA_25.1
```

For mouse, add a cached, release-pinned BioMart ortholog table and the filtered
GRCm39 BED12:

```bash
  --ortholog-table /refs/biomart/human_mouse_orthologs.tsv \
  --mouse-gene-model /refs/metagene/GRCm39/gencode_MNN/protein_coding.filtered.bed12 \
  --ensembl-release 116
```

The default retains high-confidence one-to-one orthologs. The resulting mouse
set is explicitly labelled as ortholog-derived, not as direct mouse expression
evidence.

## Standalone shared interface

Other workflows can create the documented track manifest and call:

```bash
bash common/metagene/run_metagene.sh \
  --track-manifest tracks.tsv \
  --gene-set-manifest gene_sets.tsv \
  --output-dir results/metagene \
  --gene-sets protein_coding,broadly_expressed \
  --modes tss,tes,gene_body \
  --parallel-jobs 2 \
  --threads-per-job 4
```

Track-manifest columns are:

```text
sample_id assay genome bigwig normalization normalization_detail blacklist chrom_sizes
```

Gene-set-manifest columns are:

```text
gene_set_id genome bed12 label source_release sha256 n_genes
```

Both are tab-separated. IDs must use letters, numbers, dots, underscores, or
hyphens. Gene-set BED12 checksums are verified when present.

## Outputs

Results are stored under:

```text
06_qc/metagene/<genome>/<gene_set>/<mode>/<sample>/
```

Every task writes its matrix, profile values, sorted regions, PNG/PDF profile,
PNG/PDF heatmap, and metadata. `artifacts.tsv` inventories these files, while
`multiqc_data/metagene_summary_mqc.tsv` provides compact MultiQC custom content.

Zero-only regions remain included by default. Missing coverage is interpreted
as zero only when `METAGENE_MISSING_DATA_POLICY=zero`; use `na` for tracks whose
missing intervals do not mean zero signal.

# Reference preparation and provenance

[Documentation index](README.md) | [Server installation](08_server_installation.md)

This page consolidates every reference family used by `cutnrun2tracks` 0.3.1.
It distinguishes the documented shared-server deployment from a reproducible
new build. Never infer provenance from a convenient filename such as
`annotation.gtf`; retain the upstream identity and checksum.

## Reference contract

All files for one genome must use the same assembly and chromosome vocabulary.
For the documented deployment that means UCSC-style `chr` names throughout
hg38 and mm39. The FASTA is the root of the contract:

- build the Bowtie2 index from that exact FASTA;
- derive or verify chromosome sizes against its `.fai`;
- select canonical contigs from those names;
- use assembly-matched blacklist, GTF, TSS, and cCRE intervals; and
- retain checksums before a workflow release or analysis uses them.

Do not place large references inside an immutable workflow release. References
have their own lifecycle and provenance.

## Complete configuration inventory

For each host genome used by a samplesheet, config needs:

```text
INDEX_<GENOME>
FASTA_<GENOME>
CHROM_SIZES_<GENOME>
CANONICAL_CONTIGS_<GENOME>
GTF_<GENOME>
BLACKLIST_<GENOME>
EFFECTIVE_GENOME_SIZE_<GENOME>
TSS_BED_<GENOME>          optional; derived from GTF when empty
CCRE_BED_<GENOME>         required when RUN_CCRE_ANNOTATION=true
```

Additional families are required only when enabled:

```text
SPIKEIN_INDEX, SPIKEIN_FASTA, SPIKEIN_CHROM_SIZES,
SPIKEIN_ALLOWED_CONTIGS, optional SPIKEIN_BLACKLIST

METAGENE_GENE_SET_MANIFEST and every BED12 referenced by it
```

## Documented shared-server paths

| Family | hg38 | mm39 |
|---|---|---|
| FASTA | `/opt/bioinformatics/references/hg38/hg38.fa` | `/opt/bioinformatics/references/mm39/mm39.fa` |
| Bowtie2 prefix | `/opt/bioinformatics/references/hg38/bowtie2/hg38` | `/opt/bioinformatics/references/mm39/bowtie2/mm39` |
| chromosome sizes | `/opt/bioinformatics/references/hg38/hg38.chrom.sizes` | `/opt/bioinformatics/references/mm39/mm39.chrom.sizes` |
| canonical contigs | `/opt/bioinformatics/references/cutnrun2tracks/0.2.0/hg38/hg38.canonical_contigs.txt` | `/opt/bioinformatics/references/cutnrun2tracks/0.2.0/mm39/mm39.canonical_contigs.txt` |
| GTF | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf` |
| blacklist | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.blacklist.bed` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.blacklist.bed` |
| cCRE | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed` |
| effective genome size | `2913022398` | `2654621783` |

The run-specific source of truth is the resolved config and
`00_metadata/reference_manifest.tsv`. The paths above document the established
server layout; their contents must still be checked against the shared
construction manifest.

## Host FASTA, chromosome sizes, and Bowtie2 index

### Upstream assembly choices

The reproducible UCSC sources are:

| Genome | Assembly | FASTA | chromosome sizes |
|---|---|---|---|
| hg38 | GRCh38/hg38 | [`hg38.fa.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | [`hg38.chrom.sizes`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) |
| mm39 | GRCm39/mm39, GCA_000001635.9 | [`mm39.fa.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz) | [`mm39.chrom.sizes`](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes) |

UCSC explains the hg38 sequence variants and chromosome naming in its
[hg38 bigZips README](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
and identifies the GRCm39 assembly in the
[mm39 directory](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/).
Choose one precise FASTA variant and record it; do not mix the initial hg38
FASTA with a patch/analysis-set chromosome-size file accidentally.

### Reproducible build pattern

For a new installation, adapt paths but keep the validation sequence:

```bash
GENOME=hg38
REFROOT=/opt/bioinformatics/references/$GENOME
mkdir -p "$REFROOT/bowtie2"

curl --fail --location --retry 3 \
  --output "$REFROOT/$GENOME.fa.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/$GENOME/bigZips/$GENOME.fa.gz"

curl --fail --location --retry 3 \
  --output "$REFROOT/$GENOME.upstream.chrom.sizes" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/$GENOME/bigZips/$GENOME.chrom.sizes"

gzip -t "$REFROOT/$GENOME.fa.gz"
gzip -cd "$REFROOT/$GENOME.fa.gz" > "$REFROOT/$GENOME.fa"
samtools faidx "$REFROOT/$GENOME.fa"
cut -f 1,2 "$REFROOT/$GENOME.fa.fai" > "$REFROOT/$GENOME.chrom.sizes"

diff -u "$REFROOT/$GENOME.upstream.chrom.sizes" \
        "$REFROOT/$GENOME.chrom.sizes"

bowtie2-build --threads 16 \
  "$REFROOT/$GENOME.fa" "$REFROOT/bowtie2/$GENOME"
bowtie2-inspect -n "$REFROOT/bowtie2/$GENOME" > "$REFROOT/$GENOME.index.contigs.txt"

sha256sum "$REFROOT/$GENOME.fa" "$REFROOT/$GENOME.chrom.sizes" \
  "$REFROOT"/bowtie2/$GENOME.*.bt2* > "$REFROOT/reference.sha256"
```

The `diff` must be empty. If it is not, resolve the assembly/FASTA choice rather
than editing one file until it passes. Run the equivalent block with
`GENOME=mm39`.

## Canonical-contig lists

When `CANONICAL_CHROMS_ONLY=true`, the workflow retains autosomes and X/Y from
the configured one-column list before blacklist filtering. For UCSC-style host
references:

```bash
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1}' \
  /opt/bioinformatics/references/hg38/hg38.chrom.sizes \
  > /path/to/hg38.canonical_contigs.txt

awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y)$/ {print $1}' \
  /opt/bioinformatics/references/mm39/mm39.chrom.sizes \
  > /path/to/mm39.canonical_contigs.txt

test "$(wc -l < /path/to/hg38.canonical_contigs.txt)" -eq 24
test "$(wc -l < /path/to/mm39.canonical_contigs.txt)" -eq 21
```

These lists exclude mitochondria by construction. `REMOVE_MITO` is a separate
filter for custom lists that include it. Record the chosen list and checksum.

## Blacklists

| Genome | Upstream source | Installed alias |
|---|---|---|
| hg38 | ENCODE accession [`ENCFF356LFX`](https://www.encodeproject.org/files/ENCFF356LFX/) | `hg38.blacklist.bed` |
| mm39 | Boyle Lab [`mm39-blacklist.v2.bed.gz`](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm39-blacklist.v2.bed.gz) | `mm39.blacklist.bed` |

Download, retain the compressed-source checksum, decompress without coordinate
changes, validate at least three BED columns, coordinate-sort in chromosome-size
order if necessary, and record the installed SHA-256 and interval count.

Blacklists are applied to host BAMs before coverage, peak calling, consensus,
QC, differential counting, and annotation. The scientific basis is described by
Amemiya, Kundaje, and Boyle
([Scientific Reports 2019](https://doi.org/10.1038/s41598-019-45839-z)); modern
assembly exclusions are also distributed by
[`excluderanges`](https://doi.org/10.1093/bioinformatics/btad198). Exact
filtering behavior is in
[Reference filtering and blacklists](10_references_blacklist_and_filtering.md).

## GTF, derived TSS, and cCRE references

The established annotation sources are:

| Genome | GTF | Regulatory reference |
|---|---|---|
| hg38 | GENCODE v42 `gencode.v42.primary_assembly.annotation.gtf.gz` | native GRCh38 ENCODE4 `Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz` |
| mm39 | GENCODE vM31 `gencode.vM31.primary_assembly.annotation.gtf.gz` | ENCODE3/UCSC mm10 `encodeCcreCombined.bb`, converted and lifted with `mm10ToMm39` |

[Genomic annotation](13_genomic_annotation.md) gives exact URLs, construction
steps, expected source counts, liftOver options, class vocabularies, validation
commands, and literature references. Do not duplicate or rename those sources
without retaining their provenance.

`TSS_BED_<GENOME>` is optional. When it is empty, the QC stage derives a
strand-aware one-base TSS BED from `gene` rows in the configured GTF and writes
it under:

```text
06_qc/tss_signal_profile/reference/<genome>.tss.bed
```

Supplying a shared TSS BED can save repeated derivation, but it must be generated
from the same GTF release and use the same chromosome names. Record its builder,
GTF checksum, record count, and SHA-256.

## Effective genome sizes

`EFFECTIVE_GENOME_SIZE_<GENOME>` is a positive integer used by peak callers; it
is not inferred from FASTA byte size. The documented defaults are:

```text
EFFECTIVE_GENOME_SIZE_HG38=2913022398
EFFECTIVE_GENOME_SIZE_MM39=2654621783
```

Record the value in the resolved config and methods. Changing it can alter peak
calling and requires rerunning from `peakcalling`. A custom genome requires an
explicit, scientifically justified value.

## Spike-in and competitive references

The optional spike-in branch aligns once to a composite host-plus-spike Bowtie2
index, then splits host and spike contigs. Spike names must be namespaced to
avoid collisions. Build a new composite with:

```bash
bash utilities/prepare_composite_reference.sh \
  /refs/host.fa /refs/spike.fa dm6 /refs/composite/hg38_dm6 16
```

The utility:

1. prefixes every spike contig with `<SPIKEIN_REFERENCE_ID>__`;
2. rejects host/spike name collisions;
3. creates a composite FASTA and `.fai`;
4. builds the Bowtie2 index;
5. writes the namespaced allowed-contig list; and
6. records SHA-256 values in `reference_manifest.tsv`.

For the documented dm6 deployment, existing composite prefixes are:

```text
/opt/bioinformatics/references/composite/hg38_dm6/bowtie2/hg38_dm6
/opt/bioinformatics/references/composite/mm39_dm6/bowtie2/mm39_dm6
```

and derived namespaced dm6 coordinate files live under:

```text
/opt/bioinformatics/references/cutnrun2tracks/0.2.0/dm6_namespaced/
```

`SPIKEIN_CHROM_SIZES` and `SPIKEIN_ALLOWED_CONTIGS` must contain namespaced
contigs such as `dm6__chr2L`, not original `chr2L`. The current 0.3.1 workflow
validates and records `SPIKEIN_BLACKLIST` but does not apply it to the spike BAM;
do not claim spike-reference blacklist filtering in methods.

Spike metadata (`spikein_to_host_ratio`, stage, and lot) belongs in the
samplesheet because it describes experimental libraries, not the reference.

## Metagene gene-set references

The metagene stage consumes a TSV manifest whose rows contain:

```text
gene_set_id genome bed12 label source_release sha256 n_genes
```

Build the protein-coding BED12 with
`common/metagene/prepare_metagene_reference.py`. It selects one transcript per
gene deterministically, applies canonical-contig/blacklist/length policies, and
writes metadata, exclusions, and a checksummed manifest. Build optional Human
Protein Atlas subsets with `build_hpa_reference_subset.py` from a release-pinned
[HPA download](https://www.proteinatlas.org/about/download). Mouse subsets are
ortholog-derived and require a cached, release-pinned Ensembl BioMart mapping;
the established example records HPA 25.1 and Ensembl 116.

Complete commands and semantics are in
[Metagene aggregate-signal module](06_metagene.md). Each BED12 checksum in the
gene-set manifest is revalidated before analysis.

## Provenance hierarchy

Retain all three levels:

1. **Upstream receipt:** source URL/accession, upstream filename, release,
   assembly, download date, upstream checksum, licence/terms, and publication.
2. **Construction manifest:** commands, tool versions, parameters, source
   checksums, output counts, excluded/unmapped records, output SHA-256, and
   builder date.
3. **Per-run manifest:** `00_metadata/reference_manifest.tsv`, which records
   each resolved immutable path, byte size, and SHA-256 actually used.

The run also retains `resolved_config.tsv`, `software_versions.tsv`,
`sample_manifest.tsv`, and, for composite or metagene references, their nested
manifests. A path or filename without these records is not sufficient
provenance.

## Read-only deployment audit

Before a new analysis:

```bash
cutnrun2tracks --config /absolute/path/to/config.conf --preflight-only

column -t -s $'\t' \
  /absolute/path/to/output/00_metadata/reference_manifest.tsv

bowtie2-inspect -n /opt/bioinformatics/references/hg38/bowtie2/hg38 \
  > /tmp/hg38.index.contigs.txt

diff -u \
  <(cut -f1 /opt/bioinformatics/references/hg38/hg38.chrom.sizes) \
  /tmp/hg38.index.contigs.txt
```

For large shared files, compare recorded SHA-256 values with the administrator's
reference manifest rather than trusting modification times. Stop if assemblies,
contig names, GTF/cCRE releases, or checksums differ from the intended analysis.

## Minimum acceptance checklist

- FASTA passes `samtools faidx`; Bowtie2 inspection succeeds.
- FASTA `.fai`, chromosome sizes, index contigs, and interval references share
  compatible names and lengths.
- canonical-contig counts and contents are deliberate.
- blacklist, GTF, cCRE, and optional TSS BED match the host assembly.
- compressed references pass `gzip -t`; BED/GTF records and coordinates are
  valid and non-empty.
- source and installed checksums, record counts, construction commands, and
  tool versions are retained.
- composite spike contigs are namespaced and collision-free.
- metagene manifest checksums and gene counts resolve.
- the workflow's run-specific reference manifest matches the shared manifest.

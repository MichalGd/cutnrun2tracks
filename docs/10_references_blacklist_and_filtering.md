# Reference filtering and blacklists

[Documentation index](README.md) | [Peak calling](09_peak_calling.md)

Blacklist filtering is applied to aligned observations before coverage, peak
calling, consensus, QC, differential counting, and annotation. Consequently,
blacklisted signal cannot re-enter a later host-genome result.

## Exact host blacklist files

The repository template uses assembly-specific placeholders:

```text
BLACKLIST_HG38=/refs/hg38/hg38.blacklist.bed
BLACKLIST_MM39=/refs/mm39/mm39.blacklist.bed
```

The validated shared-server deployment documented for this project uses:

| Genome | Exact configured blacklist |
|---|---|
| `hg38` | `/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.blacklist.bed` |
| `mm39` | `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.blacklist.bed` |

Those are the exact filenames, not a guarantee that every future installation
has identical contents. For a particular run, the authoritative path and
SHA-256 checksum are in `00_metadata/reference_manifest.tsv`, and the effective
key/value is in `00_metadata/resolved_config.tsv`. Blacklist assembly and
contig naming must match the Bowtie2 index, FASTA, chromosome sizes, canonical
contigs, GTF, and cCRE reference. Preflight rejects missing or empty files but
cannot infer biological provenance from a filename.

For optional dm6 spike-in on the documented server, the source file is:

```text
/opt/bioinformatics/references/dm6/dm6-blacklist.v2.bed
```

and the namespaced derivative for a composite index is:

```text
/opt/bioinformatics/references/cutnrun2tracks/0.2.0/dm6_namespaced/dm6.blacklist.bed
```

Important current limitation: `SPIKEIN_BLACKLIST` is validated and recorded,
but v0.3.1 `spikein_batch.sh` does **not** apply it to the spike BAM. Host
blacklists are applied as described below; spike-reference blacklisting should
not be claimed until the executable implementation is changed and tested.

## Why blacklist

Genome blacklists identify regions that repeatedly produce anomalous high
signal or poor mappability across many sequencing assays. Such signal can be
driven by repeats, assembly artifacts, collapsed paralogy, or other technical
effects rather than target-specific occupancy. The ENCODE blacklist analysis by
Amemiya, Kundaje, and Boyle is the primary reference
([Scientific Reports 2019, 10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z)).
The maintained [`excluderanges`](https://doi.org/10.1093/bioinformatics/btad198)
resource provides blacklist/exclusion ranges across additional assemblies,
including modern mouse assemblies.

Blacklist provenance must be checked in the site's reference manifest. Do not
silently reuse an hg19, hg38, mm10, or mm39 BED on a different assembly, even
if chromosome names appear similar.

## Filtering order and mechanism

For each host marked BAM, the workflow constructs four branches independently.
The relevant operations are:

1. retain primary mapped records and reject secondary, supplementary, QC-fail,
   and unmapped records;
2. for PE libraries, require a proper pair;
3. enforce branch-specific MAPQ (`PERMISSIVE_MIN_MAPQ=0`,
   `INTERMEDIATE_MIN_MAPQ=0`, or `MIN_MAPQ=30`);
4. retain or reject records carrying the duplicate flag according to the
   branch;
5. retain only contigs in `CANONICAL_CONTIGS_<GENOME>` when
   `CANONICAL_CHROMS_ONLY=true`, with optional mitochondrial removal controlled
   separately by `REMOVE_MITO`;
6. run `bedtools intersect -v -abam <canonical.bam> -b <blacklist.bed>` so only
   observations with no blacklist overlap survive; and
7. for PE data, name-sort, run `samtools fixmate -r`, re-require proper pairs,
   coordinate-sort, and index the result.

The PE repair step matters: if either mate is removed by blacklist or another
filter, the incomplete pair is not retained as a proper fragment. Thus the
analysis operates on complete surviving fragments rather than orphan records.
BED/BAM comparison is implemented with BEDTools; see Quinlan and Hall
([Bioinformatics 2010, 10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)).

## Duplicate handling

Picard `MarkDuplicates` marks duplicates once before the branches are made.
The workflow then provides both duplicate-retained and duplicate-removed
branches. It does not rerun an independent duplicate detector for each branch.
The target/control defaults are:

```text
TARGET_DEFAULT_DUPLICATE_POLICY=retain
CONTROL_DEFAULT_DUPLICATE_POLICY=remove
```

The samplesheet can override the policy per biological library. Retention is
often useful for sparse low-input CUT data, but high duplicate concentration
may also indicate low complexity; use the retained branch for complexity
measurement and compare retained/removed tracks as a sensitivity analysis.

## Resulting BAMs and audit files

The four branches are stored below `03_alignment/filtered/`; symlinks below
`03_alignment/analysis/` select the q30 branch used for that sample. BAMs are
indexed and content-checked before checkpoint completion. For every library,
`03_alignment/metrics/<sample>.filter_counts.tsv` records the attrition across
branches, and Picard duplicate metrics are retained separately.

When reporting results, state:

- exact blacklist path, checksum, source, version, and genome assembly;
- canonical-contig definition and whether mitochondria were removed;
- MAPQ and duplicate policy of the analysis BAM;
- numbers retained at each filter branch; and
- any change from the defaults above.

Do not describe a peak caller as performing the host blacklist removal: it
receives an already filtered analysis BAM.

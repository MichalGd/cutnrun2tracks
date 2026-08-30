# Genomic tracks and normalization

[Documentation index](README.md) | [Quality control](11_quality_control.md)

`cutnrun2tracks` produces multiple track families because library-depth,
cohort-relative, filter-sensitivity, and external-reference normalization answer
different questions. A track's suffix and directory are part of its meaning;
do not relabel one family as another.

## Coverage signal and format

All formal coverage tracks use deepTools `bamCoverage` at
`TRACK_BIN_SIZE=10` by default. For PE libraries, `--samFlagInclude 66
--extendReads` represents each retained proper fragment once from its observed
outer extent. For SE, the retained read is used without inferred fragment
extension (`SE_SIGNAL_MODE=read`). `BAMCOVERAGE_COMMON_ARGS` can add reviewed
common arguments.

Tracks are chromosome-sorted and validated before conversion:

- `.bw` is compressed/indexed and preferred for IGV, UCSC, and routine use;
- `.bedGraph.gz` is the same binned signal in compressed text form.

Retention is controlled independently by `GENERATE_COVERAGE_BIGWIGS` and
`GENERATE_COVERAGE_BEDGRAPHS`. `TRACK_STANDARD_CHROMS_ONLY=true` is consistent
with the canonical-contig filtered BAMs.

## Track-family matrix

| Family | Source BAM / observations | Scale | Meaning |
|---|---|---|---|
| analysis CPM | q30 analysis BAM; sample duplicate policy | `10^6 / L_i` | routine within-library depth-normalized signal |
| DESeq2 consensus | analysis BAM | `1 / s_i` | cohort-relative signal using consensus-count size factor |
| robust CPM permissive | MAPQ 0, duplicates retained | `10^6 / (s_i G)` | mapping/duplicate sensitivity with robust per-million units |
| robust CPM intermediate | MAPQ 0, duplicates removed | `10^6 / (s_i G)` | ambiguity sensitivity after deduplication |
| robust CPM stringent | MAPQ 30, duplicates removed | `10^6 / (s_i G)` | conservative cohort-normalized sensitivity |
| spike-in host | host analysis BAM | `K R_i / S_i` | external-reference-calibrated host signal |
| spike-control CPM | filtered spike BAM | `10^6 / S_i` | spatial/QC view of the spike reference |

Here `L_i` is the retained analysis observation count, `s_i` is a DESeq2 size
factor, `G` is defined below, `K=SPIKEIN_SCALE_TARGET`, `R_i` is the declared
spike-to-host ratio, and `S_i` is the retained spike observation count.

### A note about “raw” tracks

Version 0.3.1 does **not** publish a formal unscaled/raw bigWig family. The raw
quantitative sources are the retained BAMs and the integer consensus count
matrices. Temporary unscaled bedGraphs are implementation intermediates;
`KEEP_RAW_BEDGRAPH` can retain sorted intermediates, but that does not make them
a named normalization family. RPGC/RPKM and MACS fold-enrichment bigWigs are
also not implemented as formal `cutnrun2tracks` track families. Do not describe
CPM or DESeq2-scaled files as raw signal.

## Analysis CPM

For sample `i`, genomic-bin coverage `C_i(x)`, and retained analysis observation
count `L_i`:

```text
CPM_i(x) = C_i(x) * 1,000,000 / L_i
```

The output is `04_tracks/cpm/<sample>.CPM.{bw,bedGraph.gz}`, with the exact
count, signal unit, scale, and formula in the adjacent normalization metadata.
Each sample's samplesheet duplicate policy determines whether its source is the
q30 duplicate-retained or q30 duplicate-removed BAM.

CPM is useful for individual-library visualization and approximate comparison
when global occupancy and composition are similar. It does not correct antibody
efficiency, IP yield, batch, or a true genome-wide biological shift.

## DESeq2 consensus tracks

For each cohort, the workflow counts every target analysis BAM over the same
primary biological-support consensus. Counting uses
`GenomicAlignments::summarizeOverlaps` in `Union` mode, ignoring strand and
disallowing one alignment from being assigned to overlapping features
(`inter.feature=TRUE`). DESeq2 estimates `poscounts` size factors from the raw
integer matrix; this variant can handle sparse matrices with zeros more safely
than an all-positive geometric mean.

The track is:

```text
DESeq2Consensus_i(x) = C_i(x) / s_i
```

and is written below `04_tracks/deseq2_consensus/<cohort>/`. It is relative to
the fitted cohort and has no absolute “per million” unit. Raw counts, normalized
counts, factors, count sums, and session information are retained beside the
tracks. DESeq2 methodology is described by Love, Huber, and Anders
([Genome Biology 2014, 10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)).

## DESeq2 robust-CPM sensitivity tracks

The same cohort consensus is counted independently from each permissive,
intermediate, and stringent BAM family. For one policy:

```text
T_i = sum of sample i's consensus-region counts
G   = exp(mean(log(T_i))) across the cohort
E_i = s_i * G
robust_CPM_i(x) = C_i(x) * 1,000,000 / E_i
```

This is checked against `DESeq2::fpm(..., robust=TRUE)`. Outputs are below:

```text
04_tracks/deseq2_robust_cpm/permissive/<cohort>/
04_tracks/deseq2_robust_cpm/intermediate/<cohort>/
04_tracks/deseq2_robust_cpm/stringent/<cohort>/
```

Because each policy gets its own count matrix and size factors, a stringent
track is not created by filtering a normalized permissive track. Compare the
families to determine whether ambiguous mappings or duplicates drive a feature.
They remain cohort-relative and are not external calibration.

If a cohort lacks a usable consensus, has fewer than two targets, or contains
invalid/zero count totals, only that family is skipped by default and its reason
is recorded in `04_tracks/normalized_track_family_status.tsv`.
`REQUIRE_ALL_ENABLED_TRACKS=true` instead makes any such failure fatal.

## Spike-in normalization: end-to-end

Spike-in is disabled by default (`SPIKEIN_MODE=none`). When enabled, the
samplesheet must provide a positive `spikein_to_host_ratio`, a biologically
meaningful addition stage, and compatible lot/batch metadata. The configuration
must provide a competitive host-plus-spike Bowtie2 index, spike chromosome
sizes, allowed contigs, and reference identity.

### 1. Competitive alignment and separation

Reads are aligned once to the composite reference. Host and spike alignments
are split by their declared contigs, minimizing arbitrary assignment from two
independent alignments. The host branch proceeds through normal host filtering.
The spike branch is duplicate-marked, filtered at
`SPIKEIN_MIN_MAPQ=30`, restricted to primary/proper observations, and uses
`SPIKEIN_DUPLICATE_POLICY=remove` by default.

As documented in [Reference filtering](10_references_blacklist_and_filtering.md),
the current executable records `SPIKEIN_BLACKLIST` but does not apply it to the
spike BAM. `SPIKEIN_ALLOWED_CONTIGS` is used when separating the composite
alignment. This limitation must be disclosed for calibrated analyses.

### 2. Counts and scale factors

For sample `i`:

```text
H_i = retained host observations in the sample's analysis BAM
S_i = retained observations in the filtered spike BAM
R_i = declared spikein_to_host_ratio
K   = SPIKEIN_SCALE_TARGET (default 1,000,000)

spike_fraction_i = S_i / (H_i + S_i)
host_scale_i     = K * R_i / S_i
spike_CPM_scale_i = 1,000,000 / S_i
```

The host analysis coverage is multiplied by `host_scale_i` and written below
`04_tracks/spikein/<cohort>/`. The optional spike-reference coverage is
multiplied by `spike_CPM_scale_i` and written below
`04_tracks/spikein_control/`. Spike calibration is not applied on top of CPM or
DESeq2 scaling; it is a separate track family.

The declared ratio makes the computational factor sensitive to how much spike
material was actually introduced. A ratio of 1 for all libraries reduces the
host factor to inverse spike count times `K`. The method assumes that the spike
was added upstream of the technical variation it is intended to correct and
behaves consistently across the comparison. `spikein_stage=post_library` is
rejected as not being an upstream calibrator.

### 3. Calibration QC and failure policy

`06_qc/spikein/spikein_scaling.tsv` records the inputs, factors, spike fraction,
run median host scale, status, failures, and warnings. Defaults are:

- fail when host observations are zero;
- fail below `SPIKEIN_MIN_OBSERVATIONS_FAIL=1000` spike observations;
- warn below `SPIKEIN_MIN_OBSERVATIONS_WARN=10000`;
- warn below fraction 0.001 or above 0.20; and
- warn when the host scale is more than tenfold above or below the median of
  all valid scales in the run.

The column currently named `cohort_median_host_scale` is calculated across all
valid run rows, not separately per cohort. Use a run/samplesheet containing one
compatible spike protocol and lot context. With `ALLOW_FAILED_SPIKEIN=false`,
hard calibration failure stops the stage; when true, failed samples do not
receive spike tracks.

These checks demonstrate computational count sufficiency, not correct wet-lab
addition. Review spike lot, addition stage, cell/nuclei amount, cross-mapping,
and whether the biology could affect the calibrator. CUT&RUN calibration context
is discussed by Skene and Henikoff
([eLife 2017, 10.7554/eLife.21856](https://doi.org/10.7554/eLife.21856)) and
Meers et al. ([eLife 2019, 10.7554/eLife.46314](https://doi.org/10.7554/eLife.46314));
external-reference normalization principles are also described in ChIP-Rx
([Cell Reports 2014, 10.1016/j.celrep.2014.10.018](https://doi.org/10.1016/j.celrep.2014.10.018)).

## Choosing a family

| Goal | Starting family |
|---|---|
| inspect one library | analysis CPM |
| compare relative signal within a coherent cohort | DESeq2 consensus |
| test mapping/duplicate sensitivity | compare all three robust-CPM policies |
| compare a plausible global occupancy shift with validated spike addition | spike-in host |
| inspect spike spatial distribution and mapping artifacts | spike-control CPM |
| statistical differential binding | raw consensus counts, never bigWig values |

Always record the exact track family, BAM policy, bin size, scale-factor table,
and workflow version in figures and downstream analyses.

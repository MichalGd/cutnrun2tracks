# Peak calling and consensus construction

[Documentation index](README.md) | [Reference filtering](10_references_blacklist_and_filtering.md)

Peak calling is performed separately for every target biological library from
its final `03_alignment/analysis/<sample>.host.analysis.bam`. These BAMs have
already passed primary/proper-pair, canonical-contig, MAPQ, duplicate-policy,
and blacklist filtering. A matched control BAM is supplied whenever required.
Controls are background references, not biological target replicates.

The public defaults are:

```text
PEAK_CALLERS=seacr,macs3
PRIMARY_PEAK_CALLER=auto
REQUIRE_MATCHED_CONTROL=true
ALLOW_SHARED_CONTROLS=false
ALLOW_CONTROL_FREE_PEAKCALL=false
PEAKCALL_FAILURE_POLICY=continue
ALLOW_EMPTY_PEAKS=false
```

`target_class` in the samplesheet determines whether `narrow`, `broad`, or
both (`mixed`) peak classes are requested. Full caller outputs are retained
under `05_peaks/per_sample/`; standardized three-column BEDs are used by later
workflow stages.

## MACS3

[MACS3](https://macs3-project.github.io/MACS/docs/callpeak.html) models local
background and identifies statistically enriched treatment signal relative to
the matched control. In PE mode, `-f BAMPE` makes MACS3 use the observed paired
fragment rather than predicting and extending a single-end fragment. In SE
mode, this workflow deliberately disables model building and uses an explicit
assay-profile shift and extension.

Implemented defaults:

| Setting | PE | CUT&RUN SE | CUT&Tag SE |
|---|---|---|---|
| format | `BAMPE` | `BAM` | `BAM` |
| fragment model | observed pairs | `--nomodel --shift 0 --extsize 150` | `--nomodel --shift -75 --extsize 150` |
| narrow threshold | `-q 0.01` | same | same |
| broad threshold | `--broad --broad-cutoff 0.1` | same | same |
| duplicate handling | `--keep-dup all` | same | same |
| summits | `--call-summits` for narrow peaks | same | same |

`--keep-dup all` does not override upstream filtering: it prevents MACS3 from
applying a second duplicate rule to the already selected analysis BAM.
`MACS3_GENOME_SIZE=auto` resolves to the sum of configured chromosome sizes;
an explicit value may be supplied instead. Narrow targets produce narrowPeak
and summit outputs. Broad or mixed targets additionally use MACS3's two-level
broad procedure, in which `MACS3_QVALUE` defines stronger enriched regions and
`MACS3_BROAD_CUTOFF` joins weaker broad-domain signal.

Adjust `MACS3_QVALUE`, `MACS3_BROAD_CUTOFF`, `MACS3_CALL_SUMMITS`, or the two
SE shift/extension pairs in `config.conf`. `BOWTIE2_EXTRA_ARGS` and
`BAMCOVERAGE_COMMON_ARGS` do not change MACS3. Additional MACS3 command-line
arguments do not currently have a generic passthrough key, so changes outside
the exposed settings require a reviewed workflow change.

MACS was introduced by Zhang et al. ([Genome Biology 2008,
10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)); current
behavior and formats are defined by the linked MACS3 documentation.

## SEACR

[SEACR](https://doi.org/10.1186/s13072-019-0287-4) was designed for sparse,
low-background CUT&RUN data. It is model-free: contiguous nonzero signal blocks
are ranked by area under the curve and an empirical cutoff is selected from a
control distribution or numeric threshold.

The workflow supports SEACR for PE libraries only. It converts proper pairs to
outer-fragment coordinates, excludes fragments longer than
`SEACR_MAX_FRAGMENT=1000`, and creates a nonzero fragment-coverage bedGraph.
The target bedGraph and matched-control bedGraph are passed to SEACR with:

```text
SEACR_CONTROL_NORMALIZATION=norm
SEACR_MODE=stringent
```

`norm` asks SEACR to normalize control to target; `stringent` returns the more
selective thresholded set. When control-free calling is explicitly allowed,
`SEACR_NO_CONTROL_THRESHOLD=0.05` replaces the control bedGraph. Because the
public defaults require a control and set `SEACR_ALLOW_SE=false`, neither
fallback is used in a normal PE run.

Adjust `SEACR_MODE` (`stringent` or `relaxed`), control normalization, maximum
fragment length, or the no-control threshold in `config.conf`. Do not enable
SEACR for SE data merely to force a second caller: the implemented conversion
and validation are PE-specific.

## epic2 broad-domain sensitivity caller

[epic2](https://doi.org/10.1093/bioinformatics/btz232) bins the genome and
identifies spatially connected bins enriched over background. It is available
only when `epic2` is included in `PEAK_CALLERS` and the pinned sidecar launcher
is installed. The workflow restricts it to `broad` and `mixed` targets; it is
not run for narrow targets.

Implemented defaults are:

```text
EPIC2_BIN_SIZE=200
EPIC2_GAP_SIZE=3
EPIC2_FDR=0.05
EPIC2_FRAGMENT_SIZE=200
```

The caller receives the analysis BAM, matched control, chromosome sizes,
`--keep-duplicates`, and `--mapq 0`; the last option is intentional because
MAPQ has already been enforced upstream. PE samples add `--guess-bampe`.
Effective genome fraction is calculated as configured effective genome size
divided by the sum of chromosome sizes. The full epic2 result and a standardized
broad BED are retained.

epic2 is most useful as a broad-domain sensitivity caller for marks such as
H3K27me3, H3K9me3, or other diffuse occupancy. It should not automatically
replace MACS3: agreement and disagreement between their different spatial and
background models are informative, and caller performance remains dependent on
target biology and data quality.

## Choosing the primary caller

`PRIMARY_PEAK_CALLER=auto` resolves per cohort:

1. for broad/mixed targets, use epic2 broad peaks when epic2 is enabled;
2. for PE narrow targets, use SEACR narrow peaks when SEACR is enabled;
3. otherwise use MACS3 narrow peaks for narrow targets or MACS3 broad peaks for
   broad/mixed targets.

Set `PRIMARY_PEAK_CALLER` explicitly only when a prespecified analysis requires
one enabled, compatible caller. Non-primary callers remain useful sensitivity
outputs and are included in all-peak genomic-feature annotation.

## Per-sample failures and empty calls

Each caller/sample attempt records `SUCCESS`, `EMPTY`, or `ERROR` with its
reason and log. With `PEAKCALL_FAILURE_POLICY=continue`, a failing caller does
not terminate unrelated samples or cohorts. If the failed or empty result is
the sample's primary peak set, that sample is excluded from consensus peak
contribution but remains available for BAM-level QC and coverage. Set the
policy to `fail` for strict all-caller completion. `ALLOW_EMPTY_PEAKS` controls
whether a valid zero-interval result is acceptable.

This fault isolation prevents one low-signal sample from erasing all other
results; it does not convert the affected sample into a successful peak call.
Always review `05_peaks/peak_calling_status.tsv` and caller logs.

## Biological-support consensus

Consensus is constructed independently within each target cohort and primary
peak class. Controls never contribute peak intervals. The algorithm:

1. reads each successful target's standardized primary BED;
2. merges overlaps within that biological sample first, so one replicate cannot
   count twice at a locus;
3. sweeps all sample intervals and counts the number of distinct biological
   sample keys supporting each genomic segment; and
4. retains segments supported by at least
   `CONSENSUS_MIN_BIOLOGICAL_SAMPLES=2`.

`ALLOW_SINGLE_SAMPLE_CONSENSUS=false` prevents a one-sample cohort from being
misrepresented as replicated consensus. The resulting BED includes region IDs
and support counts and defines the fixed region universe used for FRiP,
consensus-normalized tracks, and differential binding. A successful set of
per-sample peak calls can still yield no consensus intervals when no locus
meets biological support; that is a valid, explicitly recorded skipped cohort.

## What to report

For reproducibility, report the workflow version, analysis BAM policy, caller
versions, matched-control rule, target class, all nondefault caller parameters,
per-sample caller status, primary-caller rule, consensus support threshold, and
the number of samples excluded from consensus contribution.

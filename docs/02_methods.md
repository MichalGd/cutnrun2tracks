# Methods and normalization

For PE libraries, one observation is one retained proper fragment, counted by
its first-mate record and represented over the observed outer fragment. For SE,
one observation is one retained primary read footprint; fragment size and
periodicity are not inferred.

The four filtered BAM branches are:

| Branch | MAPQ | Duplicates |
|---|---:|---|
| permissive | 0 | retained |
| intermediate | 0 | removed |
| q30 duplicate-retained | 30 | retained |
| stringent | 30 | removed |

The samplesheet duplicate policy selects one q30 branch as `analysis`.

For coverage `C_i(x)` and analysis observation count `L_i`, CPM is
`C_i(x) * 1e6 / L_i`. Consensus factors use DESeq2 `poscounts`: the relative
track is `C_i(x)/s_i`. For policy-specific peak counts with column sum `T_i`
and `G=exp(mean(log(T_i)))`, robust CPM is
`C_i(x) * 1e6/(s_i*G)`, checked against `DESeq2::fpm(robust=TRUE)`.

All factors are estimated independently for the complete target cohort key,
which includes genome, assay profile, factor, antibody ID, layout, target
class, duplicate policy, primary caller, and primary peak class.

MACS3 receives `BAMPE` for PE without shift/extension. SE receives `BAM
--nomodel` and the profile-specific provisional preset. SEACR receives nonzero
PE fragment bedGraphs and a matched control whenever required.

Optional epic2 is restricted to `broad` or `mixed` targets. It consumes the
filtered BAM and matched control, uses `--guess-bampe` for PE, configured
chromosome sizes, effective-genome fraction, bin/gap/FDR settings, and retains
its full result table beside a normalized three-column BED. epic2 is not used
for narrow targets. It is an additional sensitivity caller, not evidence that a
broad mark is automatically high quality.

## QC definitions

FastQC is run on each original technical unit, the merged biological library,
and the final trimmed reads. Library complexity is calculated from the q30
duplicate-retained BAM so duplicate removal cannot erase the multiplicity
distribution: NRF = distinct/total, PBC1 = singleton/distinct, and PBC2 =
singleton/doubleton. preseq extrapolation is descriptive and may be unavailable
for very small or low-complexity libraries. Target-control fingerprints, FRiP,
fragment-size distributions, cohort-specific Spearman matrices/PCA, and TSS
profiles provide complementary diagnostics. Optional phantompeak cross-
correlation is reported but is not promoted to a universal CUT pass/fail rule.

## Genomic-feature annotation

Every valid peak from every successful enabled caller is annotated separately;
primary consensus peaks may be included as additional entities. A peak receives
one mutually exclusive category using the configured precedence, by default:

`promoter > enhancer > exon > intron > gene_end > other_regulatory > intergenic > unclassified`.

Promoters, exons, introns, strand-aware downstream gene-end windows, nearest
genes, and signed TSS distances come from the configured GTF. Enhancer and
other-regulatory classes require the configured cCRE BED. Peaks on contigs
absent from chromosome sizes are `unclassified`. A separate all-overlaps table
preserves every feature overlap, while the exclusive assignment supports
counts, fractions, peak-covered-base fractions, and stacked horizontal plots.

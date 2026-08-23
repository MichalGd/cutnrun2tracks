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

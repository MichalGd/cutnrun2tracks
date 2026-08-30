# Limitations and pilot decisions

Before production deployment, validate per factor/antibody:

- duplicate retention versus removal;
- SEACR stringent/relaxed and MACS3 caller choice;
- CUT&RUN 0/150 and CUT&Tag -75/150 SE presets;
- the 1000-bp maximum PE insert;
- MAPQ sensitivity for repeat-associated targets;
- spike stage, count distribution, fraction, and lot consistency;
- control type and SEACR normalization mode;
- consensus support and replicate reproducibility.

The 1,000 spike-observation failure and 10,000 warning thresholds measure
numerical stability, not biological validity. Residual bacterial DNA is not
automatically a fixed external standard. Global occupancy changes require a
validated spike design or another justified normalization strategy.

Version 0.3.0 does not support SEACR for SE, mixed layouts/genomes by default,
arbitrary blocking factors in DiffBind, or multi-condition target-control
interaction models. ataqv remains disabled and, if enabled, is stored under
`experimental_ATAC_derived_ataqv` without CUT pass/fail interpretation.

epic2 is optional and isolated in a sidecar because its dependency stack is
older than the main workflow environment. Its broad domains can be sensitive
to bin size, gap allowance, effective genome fraction, control quality, and
coverage depth. MACS3 broad and epic2 results should be reviewed side by side
during pilot validation; neither is a universal default for all CUT assays.

Enhancer annotation is only as complete as `CCRE_BED_<GENOME>`. When no cCRE
reference is configured, enhancer status is explicitly `not_evaluated` and the
GTF-derived promoter/exon/intron/gene-end/intergenic classifications remain
available. Genomic categories are descriptive overlaps, not functional proof.

The CPU budget is a conservative stage-level jobs x threads check, not a cluster
scheduler or RAM estimator. It prevents obvious oversubscription when set to
`RESOURCE_CHECK_MODE=fail`, but administrators must still account for other
users, filesystem bandwidth, Java/R memory, and tool-specific transient peaks.

Consensus-normalized tracks and differential count models require every target
sample in a cohort to have nonzero counts in the retained consensus intervals.
The workflow does not silently remove zero-count samples because that would
change the intended design. In continuation mode it skips that cohort-specific
normalization/differential branch while preserving CPM tracks and other QC;
the underlying sparse or discordant data still require scientific review.

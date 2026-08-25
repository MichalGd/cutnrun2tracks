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

Version 0.2.8 does not support SEACR for SE, mixed layouts/genomes by default,
arbitrary blocking factors in DiffBind, or multi-condition target-control
interaction models. ataqv remains disabled and, if enabled, is stored under
`experimental_ATAC_derived_ataqv` without CUT pass/fail interpretation.

Consensus-normalized tracks and differential count models require every target
sample in a cohort to have nonzero counts in the retained consensus intervals.
The workflow does not silently remove zero-count samples because that would
change the intended design. In continuation mode it skips that cohort-specific
normalization/differential branch while preserving CPM tracks and other QC;
the underlying sparse or discordant data still require scientific review.

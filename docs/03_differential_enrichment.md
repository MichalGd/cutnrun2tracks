# Controls and differential enrichment

Controls have three separate uses:

1. matched background during SEACR/MACS3 peak calling;
2. target-control QC and fingerprint plots;
3. optional, separately labelled sensitivity models.

The primary analysis counts raw target fragments over the condition-unbiased
cohort consensus. Controls are not biological replicates and are not subtracted
from DESeq2 input. DESeq2Enrichment supports validated block columns and all
eligible pairwise condition contrasts.

DiffBind receives `bamReads`, `bamControl`, and `ControlID`. Its primary run
sets `bSubControl=FALSE`; the sensitivity run explicitly uses scaled control
subtraction. Version 0.1 skips DiffBind when arbitrary blocking columns are
requested, leaving the block-aware DESeq2Enrichment analysis primary.

The optional interaction model tests whether the condition change in target
signal exceeds the condition change in matched control signal. Version 0.1
requires exactly two conditions and one-to-one condition-matched controls. It
uses joint DESeq2 normalization and must be pilot-validated.

Outputs are separated into `primary_target_only`,
`sensitivity_control_subtracted`, and
`sensitivity_target_control_interaction`; only the first is primary.

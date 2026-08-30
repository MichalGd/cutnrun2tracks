# Replicates and experimental design

[Documentation index](README.md) | [Inputs and configuration](01_inputs_and_configuration.md) | [Differential binding](03_differential_enrichment.md)

This page explains how experimental units become sequencing units, biological
libraries, controls, cohorts, consensus support, and differential replicates.
Correct metadata cannot rescue a confounded experiment, but incorrect metadata
can create pseudoreplication or compare biologically incompatible libraries.

## Biological versus technical replication

**Biological replicates** are independently generated biological samples: for
example different animals, donors, independently cultured samples, or separate
experimental preparations. They capture biological variability and are the
units used for consensus support and statistical inference.

**Technical replicates** are repeated sequencing observations of the same
biological library: for example the same library distributed across lanes or
resequenced. They increase depth but do not create additional independent
biological information.

Do not use replicate numbers to describe lanes, and do not split one library
into artificial biological replicates.

## How the workflow identifies units

| Level | Samplesheet identity | Workflow handling |
|---|---|---|
| sequencing unit | one CSV row; unique `sample_id + replicate + tech_replicate` | input FASTQ pair/read |
| biological library | `sample_id + replicate` | technical FASTQs merged before trimming; one BAM/track/peak set |
| target cohort | derived scientific identity described below | consensus, normalized tracks, replicate QC, and differential analysis |

The biological output key is:

```text
<output_prefix>.bioR<replicate>
```

Biological libraries remain separate through alignment, duplicate marking,
filtering, coverage, peak calling, QC, count matrices, and differential models.
The workflow does not pool biological BAMs and then pretend the pool is a
replicate.

## Technical-replicate merging contract

Rows intended for technical merging must have the same `sample_id` and
`replicate`, distinct positive `tech_replicate` numbers, and differ only in:

```text
fastq_1, fastq_2, tech_replicate
```

Every other samplesheet field must be identical, including layout, genome,
assay profile, target, antibody, condition, control mapping, duplicate policy,
spike metadata, batch, donor, and output prefix. FASTQ paths cannot be reused in
another row.

Example: one biological target library sequenced on two lanes:

```csv
sample_id,fastq_1,fastq_2,layout,genome,assay_profile,factor,antibody_id,target_class,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_type,control_id,analysis_duplicate_policy,spikein_to_host_ratio,spikein_stage,spikein_lot,batch,donor,output_prefix
WT_H3K27ac,/data/lane1_R1.fastq.gz,/data/lane1_R2.fastq.gz,PE,hg38,cuttag,H3K27ac,AB_K27AC,broad,WT,vehicle,keratinocyte,1,1,FALSE,none,WT_IgG,retain,,,,B1,D1,WT_H3K27ac
WT_H3K27ac,/data/lane2_R1.fastq.gz,/data/lane2_R2.fastq.gz,PE,hg38,cuttag,H3K27ac,AB_K27AC,broad,WT,vehicle,keratinocyte,1,2,FALSE,none,WT_IgG,retain,,,,B1,D1,WT_H3K27ac
```

These rows become one biological library, `WT_H3K27ac.bioR1`. They count once
for consensus support and once in a differential design.

## Biological-replicate representation

Independent replicates use the same biological `sample_id` with different
positive `replicate` values, or deliberately distinct safe sample IDs with
their own output prefixes. Each target replicate normally has a corresponding
control with the same replicate number.

```csv
WT_H3K27ac,...,1,1,FALSE,none,WT_IgG,...,D1,WT_H3K27ac
WT_H3K27ac,...,2,1,FALSE,none,WT_IgG,...,D2,WT_H3K27ac
WT_IgG,...,1,1,TRUE,igg,,...,D1,WT_IgG
WT_IgG,...,2,1,TRUE,igg,,...,D2,WT_IgG
```

The abbreviated rows above illustrate identities only; production files must
contain the exact full header and every required field.

## Matched-control design

A target's `control_id` names a control row's `sample_id`. With the default
policy, exactly one compatible control must resolve for the same biological
replicate.

### Ordinary replicate-matched controls

Target and control must agree in:

- `replicate`;
- genome and PE/SE layout;
- CUT&RUN/CUT&Tag assay profile;
- condition, treatment, cell type, and batch; and
- spike-in ratio, stage, and lot when spike-in is enabled.

The target row uses `is_control=FALSE`, `control_type=none`, and
`target_class=narrow|broad|mixed`. The control row uses `is_control=TRUE`,
`target_class=control`, `control_type=igg|input|mock`, and an empty `control_id`.

Control factor and antibody identifiers describe the control itself; they need
not equal the target factor/antibody. Control duplicate policy can also differ
from the target policy.

### Shared controls

Shared controls are disabled by default. If an experimental design genuinely
uses one control across targets or replicates:

```text
ALLOW_SHARED_CONTROLS=true
```

A fallback shared control must have `condition=shared`. It still must match
genome, assay, layout, treatment, cell type, batch, and spike metadata. Every
reuse appears in `cohort_membership.tsv`.

Shared controls reduce sequencing cost but can create dependence among
target-control summaries and cannot substitute for independent controls in the
target-control interaction model. Document the experimental reason rather than
enabling sharing to bypass an unresolved mapping error.

### Control-free peak calling

`ALLOW_CONTROL_FREE_PEAKCALL` and caller-specific thresholds permit special
control-free designs, but the default requires matched controls. A control-free
run changes the peak-calling estimand, disables control-based QC, and makes
control-aware sensitivity models unavailable. It should be a deliberate design,
not a post hoc repair for missing libraries.

## Cohort construction

Targets form the same cohort only when all of these agree:

```text
genome
assay_profile
factor
antibody_id
layout
target_class
analysis_duplicate_policy
primary_peak_caller
primary_peak_class
```

When spike-in is enabled, the identity also includes spike-in mode, reference,
stage, and lot. The cohort ID has a readable prefix plus a hash of the complete
identity; inspect `cohort_manifest.tsv` rather than reverse-engineering the
filename.

Consequences:

- different factors are never normalized or tested together;
- different antibody IDs form separate cohorts even if factor names match;
- narrow, broad, and mixed targets can select different primary callers/classes;
- retain/remove duplicate policies do not silently mix; and
- PE and SE signal units do not share a cohort.

If nominal replicates unexpectedly form separate cohorts, compare their full
cohort keys. Do not erase a meaningful antibody, assay, or policy difference
merely to force them together.

## Consensus support

Consensus uses only successful primary target peak sets. Support is counted
across distinct biological sample keys after technical merging. Controls never
contribute peaks or support.

Default behavior:

```text
CONSENSUS_MIN_BIOLOGICAL_SAMPLES=2
ALLOW_SINGLE_SAMPLE_CONSENSUS=false
CONSENSUS_USE_PRIMARY_CALLER=true
```

A cohort with three successful peak samples can retain intervals supported by
at least two, even if a fourth biological sample was excluded after a caller
failure. The summary and `excluded_peak_samples.tsv` make this denominator
visible.

This support rule is not formal IDR. A peak present in two low-quality samples
is not automatically trustworthy, and a real condition-specific peak may fail
an all-sample support threshold. Review per-sample QC and caller sensitivity
sets alongside the consensus.

## Differential eligibility and models

A condition enters modeling only when it has at least:

```text
DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION=2
```

independent target biological libraries in the cohort. Technical units do not
increase this count. At least two eligible conditions are required. For a new
study, three or more biological replicates per condition is generally more
robust than the executable minimum, especially with heterogeneous primary
material, but power should be justified for the expected effect, variability,
and multiple-testing burden.

The primary direct model is:

```text
~ condition
```

or, with prespecified blocks:

```text
~ block_1 + block_2 + ... + condition
```

Set samplesheet `batch` and `donor`, then list the required fields in
`DIFFERENTIAL_BLOCK_COLUMNS`. Every modeled sample must have a non-placeholder
value for each requested block, each block must have at least two levels, and
the design matrix must be full rank.

### Common valid designs

| Experimental question | Metadata/design |
|---|---|
| independent WT versus KO | balanced biological replicates; `~ condition` |
| treatment with samples processed in multiple batches | distribute every condition across batches; `~ batch + condition` |
| paired pre/post samples from donors | every donor contributes both conditions; `~ donor + condition` |
| multi-condition experiment | replicated conditions; workflow exports every eligible pairwise contrast |

### Confounded designs that cannot be repaired statistically

- every control sample is in batch 1 and every treatment sample in batch 2;
- each donor appears in only one condition when donor and condition effects must
  be separated;
- all replicates of one condition use a different antibody lot recorded as a
  separate antibody cohort;
- two lanes from one library are entered as two biological replicates; or
- one control library is copied into multiple rows and presented as independent
  controls.

When condition is perfectly confounded with batch/donor, removing the block can
produce a model but does not recover the missing experimental information.

## Primary and control-aware analyses

The primary target-only analysis counts raw target fragments/reads over the
fixed cohort consensus. Controls are not target replicates and are not
subtracted from its DESeq2 count matrix.

Control-aware outputs are separate sensitivity analyses:

- DiffBind scaled control subtraction changes the estimand and introduces
  control sampling noise;
- the target-control interaction requires exactly two conditions, one-to-one
  non-reused condition-matched controls, and a full-rank joint design; and
- primary DiffBind can use controls for greylisting without count subtraction.

Plan control replication for the sensitivity question you intend to ask. One
shared IgG may support peak background estimation but cannot satisfy the
one-to-one interaction design.

## Batch, donor, and antibody effects

- Randomize conditions across extraction, library, sequencing, and processing
  batches where possible.
- Keep batch labels at the biological-library level. Two sequencing lanes of
  the same library cannot carry different biological batch labels if they are to
  merge.
- Record donor for paired designs before examining PCA.
- Use a stable `antibody_id` to distinguish reagents/lots whose effects should
  not be assumed exchangeable. Different IDs form different cohorts rather than
  being adjusted as a block within one cohort.
- Inspect replicate Spearman heatmaps and PCA, but do not choose or remove
  blocking variables solely to maximize significance.

The DiffBind wrapper skips arbitrary block designs in 0.3.1; the direct
DESeq2Enrichment path implements the additive blocks and remains primary.

## PE/SE, assay, target-class, and duplicate-policy design

Keep mixed genomes and mixed layouts in separate runs unless the expert
overrides are necessary and validated. PE units are physical fragments; SE
units are retained reads. They are not equivalent depths.

CUT&RUN and CUT&Tag profiles can coexist as metadata, but assay profile is part
of the cohort identity. Similarly:

- use `narrow` for punctate binding;
- use `broad` for domain-like marks;
- use `mixed` when both caller classes are scientifically relevant; and
- choose duplicate retain/remove policy before inspecting results, then use the
  alternative branch as sensitivity evidence when needed.

Target class affects primary caller resolution. If epic2 is enabled, broad or
mixed targets can select it automatically; narrow targets cannot use epic2 as
primary.

## Spike-in experimental design

Spike normalization is valid only when spike material, amount, stage, and lot
are experimentally interpretable and consistently recorded. The cohort identity
separates spike stage and lot. A numerical factor cannot correct an unrecorded
change in spike addition, species mixture, recovery, or batch.

Plan sufficient spike observations for every modeled target and control. Failed
spike QC can invalidate calibrated tracks and spike-normalized differential size
factors even when host alignment succeeds.

## Design audit before sequencing or analysis

Before a full run:

1. identify the independent biological unit and assign `replicate` accordingly;
2. identify lanes/resequencing of the same library as `tech_replicate`;
3. ensure every intended target replicate has the correct control replicate;
4. balance condition across batch, donor, spike lot, and other technical factors;
5. confirm antibody IDs, factor, target class, assay, layout, and duplicate
   policy produce the intended cohorts;
6. ensure each intended differential condition meets the replicate minimum;
7. prespecify blocks and contrasts before inspecting results; and
8. run the metadata plan:

```bash
cutnrun2tracks --config /absolute/path/to/config.conf --plan

column -t -s $'\t' <OUTPUT_DIR>/00_metadata/sample_manifest.tsv | less -S
column -t -s $'\t' <OUTPUT_DIR>/00_metadata/cohort_manifest.tsv | less -S
column -t -s $'\t' <OUTPUT_DIR>/00_metadata/cohort_membership.tsv | less -S
```

Check the reported numbers of sequencing units, biological libraries, and
target cohorts. Resolve discrepancies before starting preprocessing.

## Interpretation checklist after the run

- verify every technical unit contributed to the intended biological library;
- review target-control fingerprints and control mappings;
- inspect per-sample peak status and exclusions from consensus;
- compare sequencing depth, filtering, duplicates, NRF/PBC, FRiP, and peak
  counts without imposing universal assay thresholds;
- inspect within-cohort correlation and PCA for outliers and batch structure;
- confirm eligible conditions, design formula, rank, reference condition, and
  contrast direction;
- interpret sensitivity analyses as different estimands, not replicated tests;
  and
- document exclusions and deviations from the prespecified design.

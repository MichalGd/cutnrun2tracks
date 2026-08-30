# Differential binding analysis

[Documentation index](README.md) | [Replicate inputs](01_inputs_and_configuration.md)

Differential binding tests whether raw target fragment/read counts over a fixed
consensus-region universe differ between biological conditions. It does not
test bigWig or bedGraph values, and it does not treat controls as target
replicates. `cutnrun2tracks` separates the primary target-only estimand from two
control-aware sensitivity analyses.

## Shared region universe and count input

Each target cohort is isolated by genome, assay profile, factor, antibody ID,
layout, target class, duplicate policy, primary caller, and primary peak class.
The cohort's primary per-sample peak sets are combined into a consensus that
requires support from distinct biological samples. This condition-unbiased
consensus is fixed before testing.

The normalized-tracks stage counts every target analysis BAM over those regions
with `GenomicAlignments::summarizeOverlaps`, `mode="Union"`, ignored strand,
and `inter.feature=TRUE`. The retained non-negative integer matrix is:

```text
04_tracks/deseq2_consensus/<cohort>/tables/raw_counts.tsv.gz
```

It is generated when needed by either normalized tracks or differential
analysis. Conditions with fewer than
`DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION=2` biological libraries are excluded
from modeling; at least two eligible conditions are required. Technical FASTQ
units have already been merged and never increase replicate count. Controls do
not contribute consensus intervals or target replicate counts.

## Analysis families

| Family | Method | Role | Use of matched controls |
|---|---|---|---|
| `primary_target_only/deseq2_enrichment` | direct DESeq2 | primary | excluded from model |
| `primary_target_only/diffbind` | DiffBind with DESeq2 | primary peer | attached for greylisting; no count subtraction |
| `sensitivity_control_subtracted/diffbind` | DiffBind with DESeq2 | sensitivity | scaled subtraction enabled |
| `sensitivity_target_control_interaction/deseq2` | joint DESeq2 interaction | sensitivity | modeled as paired control libraries |

Only target-only analyses are designated primary. Agreement among methods can
support robustness, but disagreement is diagnostic: the methods have different
counting, normalization, background, and design behavior.

## Primary direct DESeq2 model

The design is:

```text
~ condition
```

or, when `DIFFERENTIAL_BLOCK_COLUMNS` is set:

```text
~ block_1 + block_2 + ... + condition
```

Typical blocks are donor, batch, or another prespecified nuisance factor in the
samplesheet. Every requested block must exist, have no `.` values in eligible
samples, and contain at least two levels. The complete design matrix must be
full rank; confounding (for example one batch unique to each condition) is an
error rather than silently dropped. The model is additive: v0.3.1 does not
expose arbitrary condition-by-block interactions in this primary path.

DESeq2 fits a negative-binomial generalized linear model per consensus region,
estimates dispersions across regions, and applies multiple-testing correction.
With the default `DIFFERENTIAL_NORMALIZATION=deseq2`, size factors are estimated
using `type="poscounts"`, which is appropriate for a sparse peak-count matrix
containing zeros. With `DIFFERENTIAL_NORMALIZATION=spikein`, the workflow sets:

```text
effective_i = spike_observations_i / spikein_to_host_ratio_i
size_factor_i = effective_i / geometric_mean(effective)
```

This requires valid spike-in QC for every modeled target and changes only the
normalization factor, not the raw count matrix or design.

Condition ordering follows first appearance among eligible samples unless
`DIFFERENTIAL_CONDITION_ORDER` is specified. An optional
`DIFFERENTIAL_REFERENCE_CONDITION` is moved to the reference position. Every
pairwise contrast among eligible conditions is exported as
`numerator_vs_reference`; positive log2 fold change means higher target signal
in the numerator.

A region is called significant only when:

```text
padj <= DIFFERENTIAL_ALPHA             (default 0.05)
abs(log2FoldChange) >= DIFFERENTIAL_MIN_ABS_LOG2FC  (default 0)
```

`padj=NA` rows remain in the all-results table but are not significant. The
workflow does not apply post-hoc fold-change shrinkage. VST PCA uses
`blind=FALSE`; the workflow also writes a dispersion plot, raw/normalized
counts, size factors, comparison summaries, serialized DESeq2 object, and R
session information. See Love et al.
([Genome Biology 2014, 10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8))
and the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## Primary DiffBind peer analysis

The DiffBind samplesheet contains each target analysis BAM, its matched control
analysis BAM, and the same fixed consensus BED. In the validated default design,
every target must therefore resolve to a readable matched control. DiffBind:

1. imports the fixed consensus rather than constructing a new peak union;
2. disables its additional genome blacklist because host BAMs and the consensus
   are already blacklist-filtered upstream;
3. enables control-derived greylists;
4. calls `dba.count(..., bUseSummarizeOverlaps=TRUE,
   bSubControl=FALSE, bScaleControl=FALSE)`; and
5. analyzes replicated condition contrasts with DiffBind's DESeq2 method.

This is still target-only testing: the control informs greylisting but its count
is not subtracted. DiffBind manages its own count/normalization object, so its
statistics need not equal the direct DESeq2 results. It applies
`DIFFERENTIAL_ALPHA` when exporting significant results; the direct model's
`DIFFERENTIAL_MIN_ABS_LOG2FC` threshold is not separately applied by the
current DiffBind script.

If arbitrary block columns are configured, v0.3.1 writes a DiffBind `SKIPPED`
record because that wrapper does not map them into a DiffBind formula. The
block-aware direct DESeq2 analysis remains primary. See the
[DiffBind Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
and its reference manual for `dba.count`, `dba.blacklist`, and DESeq2 analysis.

## Sensitivity 1: control-subtracted DiffBind

When `RUN_CONTROL_SUBTRACTED_SENSITIVITY=true`, the same DiffBind analysis is
repeated with:

```text
bSubControl=TRUE
bScaleControl=TRUE
```

This tests target signal after DiffBind's scaled matched-control subtraction.
It may reduce condition-dependent background, but also introduces control
sampling noise and changes the estimand. It requires complete matched controls
and is skipped for arbitrary block designs by the same wrapper rule. It should
not replace the primary target-only result merely because it yields more or
fewer significant regions.

## Sensitivity 2: target-control interaction

When `RUN_TARGET_CONTROL_INTERACTION=true`, the workflow jointly counts targets
and uniquely matched controls over the target consensus and fits:

```text
~ library_type + condition + library_type:condition
```

`library_type` uses control as baseline. The interaction coefficient therefore
tests:

```text
(target change between conditions) - (control change between conditions)
```

The current implementation requires exactly two eligible conditions,
one-to-one target/control matching, no reused control, condition-matched control
rows, and a full-rank design. Size factors are estimated jointly with DESeq2
`poscounts`. Significant interaction results use adjusted P value
`DIFFERENTIAL_ALPHA`; the direct model's minimum absolute fold-change threshold
is not applied here. This is a pilot sensitivity model, not the default primary
analysis.

## Control policy keys

`DIFFERENTIAL_CONTROL_MODE=peak_calling_only` documents that the primary model
does not subtract control. `DIFFERENTIAL_SUBTRACT_CONTROL=false` is a safety
guard; enabling direct blanket subtraction is rejected. Control use is instead
explicitly isolated into the named DiffBind and interaction sensitivity
families. `REQUIRE_CONDITION_MATCHED_CONTROLS=true` protects the interaction
design.

## Outputs

The root is `08_differential/<cohort>/<peak_class>/`. The direct DESeq2 path
contains:

```text
primary_target_only/deseq2_enrichment/
|-- raw_counts.tsv.gz
|-- normalized_counts.tsv.gz
|-- size_factors.tsv
|-- comparison_summary.tsv
|-- pca.tsv and pca.png
|-- dispersion.png
|-- deseq2_object.rds
|-- session_info.txt
`-- comparisons/<numerator>_vs_<reference>/
    |-- all.tsv.gz
    `-- significant.tsv.gz
```

DiffBind paths retain their generated samplesheet, complete/significant contrast
TSVs, R object, summary, and session. The interaction path retains complete and
significant gzipped tables and its DESeq2 object. The report stage indexes all
`SUCCESS`, `SKIPPED`, `FAILED`, `DISABLED`, and unavailable variants rather than
showing only successful comparisons.

The annotation stage adds adjacent annotated copies to gzipped tables that
contain `region_id`; original statistical results remain unchanged. Current
DiffBind contrast TSVs are uncompressed and are not automatically annotated.

## Interpretation and sensitivity checklist

1. Confirm the contrast direction, eligible sample count, condition order,
   design formula, rank, normalization, thresholds, and software/session info.
2. Review sequencing depth, filtering, complexity, target/control separation,
   per-sample peak status, consensus support, FRiP, within-condition correlation,
   PCA, and donor/batch balance before interpreting P values.
3. Use direct DESeq2 as the block-aware primary model; treat primary DiffBind as
   a peer implementation, not a technical replicate.
4. Compare primary results with control-subtracted and interaction results only
   when their matching assumptions are satisfied.
5. Investigate sign reversals or sensitivity-only discoveries for background,
   low counts, confounding, normalization failure, or a genuine change in
   control signal.
6. Zero significant regions is a valid statistical outcome. It is not a reason
   to relax thresholds after viewing results.
7. Nearest-gene or cCRE annotations provide context but do not establish a
   regulated gene or mechanism.

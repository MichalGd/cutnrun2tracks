# Quality control

[Documentation index](README.md) | [Tracks and normalization](12_tracks_and_normalization.md)

Quality control is divided here into **sequencing QC**, which asks whether the
FASTQs and alignments are technically usable, and **assay-performance QC**,
which asks whether the expected target-specific chromatin signal is present and
reproducible. `QC_THRESHOLDS_MODE=descriptive` is the implemented policy: the
workflow reports evidence but does not silently exclude a library because one
metric crossed a generic threshold.

## Sequencing QC

### FASTQ integrity and FastQC

When `RUN_FASTQC=true`, FastQC is run at three levels:

1. each original technical sequencing unit when
   `RUN_FASTQC_PER_TECHNICAL_UNIT=true`;
2. the concatenated raw biological-library FASTQ(s); and
3. the final trimmed biological-library FASTQ(s).

This preserves lane/unit diagnostics while also showing what the aligner
actually receives. FastQC modules include basic statistics, sequence count and
length, per-base and per-sequence quality, base composition, GC distribution,
per-base N content, sequence-length distribution, duplication, overrepresented
sequences, adapter content, and per-tile quality when the input headers support
it. FastQC's warning/failure icons assume a relatively random library and must
be interpreted in assay context; enriched libraries can legitimately violate
that assumption. See the [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
and its guidance on [interpreting module flags](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/2%20Basic%20Operations/2.2%20Evaluating%20Results.html).

Outputs are below `01_fastq_qc/raw_units/`, `01_fastq_qc/raw/`, and
`01_fastq_qc/trimmed/`. `RUN_MULTIQC=true` creates the preprocessing summary
and later the final unified MultiQC report.

### Adapter trimming

When `TRIM_ADAPTERS=true`, Trim Galore/cutadapt processes each merged library.
The defaults use paired mode for PE, single mode for SE,
`THREADS_TRIMGALORE=8`, and `MIN_TRIMMED_LENGTH=20`; additional arguments come
from `TRIMGALORE_EXTRA_ARGS`. Review the retained read/pair count, length
filtering, adapter sequence and base removal, and the before/after FastQC
modules. A large loss can indicate short inserts, adapter dimers, low-quality
cycles, or an inappropriate adapter choice.

### Alignment and filter retention

Bowtie2 logs record overall and concordant alignment summaries. The public
defaults are end-to-end, very-sensitive, PE insert range 10--1000, dovetail
allowed, mixed and discordant pairs disabled, best alignment, and seed 0.
`samtools flagstat` and `samtools stats` are retained for the analysis BAM.

Every sample also has a filtering-count table covering q0/q30 and
duplicate-retained/removed branches. It reveals attrition from MAPQ, canonical
contigs, blacklist removal, proper-pair repair, and duplicate policy. Review
these counts before interpreting a shallow final BAM; the mechanism is detailed
in [Reference filtering](10_references_blacklist_and_filtering.md).

### Read/fragment length

For PE libraries, `RUN_FRAGMENT_QC=true` reports observed template lengths from
the analysis BAM up to `FRAGMENT_PLOT_MAX_BP=1000`. The histogram can expose
adapter-sized inserts, an unexpected long-fragment tail, or library-to-library
shifts. It is descriptive: CUT&RUN and CUT&Tag profiles depend on target and
protocol, and the workflow does not require an ATAC-like nucleosomal pattern.
SE data have no inferred fragment-length QC.

## Assay-performance QC

### Duplicate metrics and library complexity

Picard duplicate metrics report marked duplicate burden. The workflow's
coordinate-complexity metrics deliberately use the **q30 duplicate-retained**
BAM, because calculating them after deduplication would erase the multiplicity
being measured:

- **NRF** = distinct genomic signal locations / total signal units;
- **PBC1** = singleton locations / distinct locations;
- **PBC2** = singleton locations / doubleton locations (`Inf` when no
  doubleton exists).

PE coordinates are outer fragments; SE coordinates are read footprints.
Picard and coordinate multiplicity are related but not identical definitions,
so their numbers need not match exactly.

When `RUN_PRESEQ=true`, `preseq lc_extrap -B` estimates the number of distinct
molecules expected with additional sequencing. It is descriptive and may fail
for a tiny or extremely sparse multiplicity distribution; the pipeline records
a warning and continues. The method is described by Daley and Smith
([Nature Methods 2013, 10.1038/nmeth.2375](https://doi.org/10.1038/nmeth.2375)).

### Cross-correlation

When `RUN_CROSS_CORRELATION=true`, phantompeakqualtools processes a tagAlign
made from the q30 duplicate-retained BAM. It reports strand cross-correlation,
including normalized strand coefficient (NSC), relative strand correlation
(RSC), estimated fragment-length peak, read-length phantom peak, and associated
plots/status. NSC and RSC measure enrichment independently of a called peak
set; the [ENCODE tool record](https://www.encodeproject.org/software/phantompeakqualtools/)
links their definitions and reference publication.

Cross-correlation thresholds developed for ChIP-seq transcription factors are
not automatically valid for sparse CUT assays or broad histone domains. This
module is off by default and should remain descriptive unless a project has
validated target-specific criteria.

### Peak-call and consensus diagnostics

`05_peaks/peak_calling_status.tsv` records every sample/caller/class attempt,
its status, peak count, and failure reason. Consensus status records total
target samples, successful primary-peak contributors, exclusions, support
threshold, and consensus region count. Review these before treating a zero-peak
or excluded sample as merely a software warning.

A library may have valid alignments and coverage but no reproducible peaks.
Conversely, one sample with many idiosyncratic peaks can fail to contribute to
a supported consensus. Both outcomes are scientifically informative.

### FRiP

For each target, FRiP is the fraction of analysis signal units overlapping its
cohort's **primary biological-support consensus**, not that sample's own peak
set. PE fragments are counted once using outer coordinates; SE reads are
counted once. Using a common region universe makes replicates comparable within
the cohort, but FRiP still depends on caller, target class, consensus support,
MAPQ, duplicate policy, and mark breadth. It is not comparable without those
definitions.

### Target-control fingerprint

For every target with a matched control, deepTools `plotFingerprint` compares
cumulative genome-wide signal distributions from the two analysis BAMs. Clear
separation supports enrichment; weak separation can reflect low target signal,
high background, poor antibody specificity, or a biologically diffuse target.
The matched control is not included as a target replicate.

### Replicate correlation and PCA

When `RUN_REPLICATE_CORRELATION=true`, target analysis BAMs in a cohort are
summarized in genomic bins with `multiBamSummary bins`. The workflow produces a
Spearman correlation heatmap and PCA. Interpret within-condition replicate
agreement together with experimental design: separation by condition can be
desired, whereas separation by lane, donor imbalance, or processing batch may
indicate confounding. Correlation across an overwhelmingly empty genome may
also obscure differences visible in consensus-region counts.

### TSS signal profile

When `RUN_TSS_SIGNAL_PROFILE=true`, the workflow makes an aggregate profile of
analysis CPM signal around strand-aware TSS coordinates, using
`TSS_BED_<GENOME>` or TSSs derived from GTF gene records. Defaults are 3 kb
upstream and downstream. This is a descriptive average profile, **not** the
ATAC-seq TSS-enrichment ratio. Many chromatin targets are not expected to peak
at promoters.

### Optional ataqv

`RUN_ATAQV_QC=false` by default. If explicitly enabled, ataqv receives target
analysis BAMs, primary peaks, the configured blacklist as an exclusion set,
and a TSS extension. Its metrics and optional viewer originate in an ATAC-seq
QC framework. The versioned ataqv JSON/viewer is the complete output and
includes its standard read-accounting, mapping/duplicate, autosomal and
mitochondrial, fragment-length, peak-overlap/FRiP, and TSS-enrichment fields
when calculable. The CUT workflow does not promote those fields into automatic
thresholds. They are experimental supplementary diagnostics, not validated
universal CUT pass/fail criteria.

### Spike-in performance

When spike-in is enabled, the QC table reports host and spike observations,
spike fraction, declared spike-to-host ratio, host and spike scaling factors,
status, reason, and warnings. Defaults fail below 1,000 spike observations,
warn below 10,000, warn outside a 0.001--0.20 spike fraction, and warn when a
host scale differs by more than tenfold from the run median. Full formulas and
experimental cautions are in [Tracks and spike-in normalization](12_tracks_and_normalization.md).

### Metagene profiles

Optional TSS, TES, and scaled gene-body profiles provide biological context for
chosen gene sets and track families. They are aggregate signal visualizations,
not independent library-quality scores. See [Metagene analysis](06_metagene.md).

## QC output map

| Question | Main evidence | Typical location |
|---|---|---|
| Are raw/trimmed reads sound? | FastQC, Trim Galore/cutadapt, preprocessing MultiQC | `01_fastq_qc/`, `logs/preprocess/` |
| Did alignment/filtering retain enough signal? | Bowtie2 log, flagstat, stats, branch counts | `03_alignment/metrics/`, `06_qc/alignment_and_complexity/` |
| Is the library complex? | Picard, NRF/PBC, preseq | same directories plus `logs/qc/` |
| Is fragment structure plausible? | PE fragment histogram | `06_qc/fragment_length_and_periodicity/` |
| Is enrichment detectable without peaks? | NSC/RSC, target-control fingerprint | `06_qc/fragment_length_and_periodicity/`, `06_qc/controls/` |
| Is signal in reproducible peaks? | caller status, consensus status, FRiP | `05_peaks/`, `06_qc/frip_and_peak_reproducibility/` |
| Do replicates agree? | Spearman heatmap, PCA | `06_qc/correlation_pca_fingerprint/` |
| Is promoter-centered signal present? | descriptive TSS profile | `06_qc/tss_signal_profile/` |
| Is calibration credible? | spike-in scaling/status table | `06_qc/spikein/` |

## Final review

The unified report is `10_reports/cutnrun2tracks_multiqc_report.html`; stable
TSV tables remain available for programmatic review. Before accepting or
excluding a library, record the reason and evaluate at least sequencing
quality, alignment/filter retention, complexity, target-control separation,
peak status, consensus support, FRiP, replicate structure, and any spike-in
warnings together. A FastQC or preseq warning alone is not a sufficient assay
decision.

The principal aggregation and coverage tools are MultiQC
([Bioinformatics 2016, 10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354))
and deepTools
([Nucleic Acids Research 2016, 10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)).
Their plots summarize the inputs supplied by this workflow; they do not change
the filtering or statistical definitions documented above.

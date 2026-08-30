# Configuration files

Copy `config.conf.template` to an untracked `config.conf`, then set every input,
output, and reference path used by the selected genome. The parser accepts only
the keys present in the template plus assembly-specific reference keys matching
`INDEX_*`, `FASTA_*`, `CHROM_SIZES_*`, `CANONICAL_CONTIGS_*`, `GTF_*`,
`BLACKLIST_*`, `EFFECTIVE_GENOME_SIZE_*`, `TSS_BED_*`, or `CCRE_BED_*`.

Sample-level properties occur only in the samplesheet: FASTQ paths, genome,
layout, assay profile, target, antibody, condition, replicate, control mapping,
duplicate policy, and optional spike metadata. Reference paths and run policy
occur only in `config.conf`; the blacklist is resolved from `BLACKLIST_<GENOME>`
and injected into the internal manifest. There is no run-wide
`ASSAY_PROFILE` setting.

`RUN_METAGENE=true` additionally requires `METAGENE_GENE_SET_MANIFEST` to name
a tab-separated gene-set manifest prepared with the utilities under
`common/metagene/`. The plotting stage consumes existing CPM or spike-in
bigWigs and does not perform normalization.

The examples contain illustrative `/data` and `/refs` paths and are not directly
executable. `cutrun_se.csv` requires `PEAK_CALLERS=macs3`.

For broad domains, add `epic2` to `PEAK_CALLERS` only after installing the
versioned sidecar from `environment.epic2.yml`. With
`PRIMARY_PEAK_CALLER=auto`, broad/mixed targets use epic2 when enabled, narrow
PE targets prefer SEACR, and remaining targets use MACS3. All enabled successful
caller peak sets are retained and annotated.

`RUN_MULTIQC=true` generates the preprocessing FastQC aggregation and the final
unified report under `10_reports/`. The final report scans retained standard-tool
logs and adds CUT-specific custom tables for observations, filtering, peak calls,
consensus, normalization, spike-in QC, metagene, and differential status.

Peak-caller fault handling is explicit. `PEAKCALL_FAILURE_POLICY=continue`
records each enabled caller as `SUCCESS`, `EMPTY`, or `ERROR`, keeps processing
unaffected samples, and excludes only failed/empty primary peak contributions
from consensus construction. A cohort continues when at least
`CONSENSUS_MIN_BIOLOGICAL_SAMPLES` successful primary peak sets remain. Use
`PEAKCALL_FAILURE_POLICY=fail` for strict fail-fast operation. In either mode,
`ALLOW_EMPTY_PEAKS` controls whether a zero-peak caller result is accepted or
flagged as a problem; it never fabricates peaks.

Normalized-track fault handling uses `REQUIRE_ALL_ENABLED_TRACKS`. Its default
`false` value records and skips only a cohort/track family whose consensus is
unavailable or whose consensus counts cannot be normalized, then continues
unaffected work. Zero-count samples are reported, never silently discarded.
Set it to `true` when every requested normalized family is mandatory.

Preprocessing parallelism is controlled at two levels:

- `QC_SAMPLE_PARALLEL_JOBS` limits concurrent biological-library workers;
- `THREADS_FASTQC` is passed to each FastQC process;
- `THREADS_TRIMGALORE` is passed as `trim_galore --cores` for each worker.

Their products approximate the maximum preprocessing CPU demand. For example,
the template's `8 x 10` FastQC and `8 x 8` Trim Galore settings suit a large
shared server but should be reduced for smaller hosts or concurrent heavy runs.

`TOTAL_CPU_BUDGET` provides a run-wide ceiling check. Preflight compares each
long stage's configured jobs x threads against it and writes
`00_metadata/resource_budget.tsv`; `RESOURCE_CHECK_MODE=fail` prevents an
overcommitted run from starting. Separate job limits are available for
alignment/filtering, tracks, peak calling, spike-in, normalized tracks,
differential analysis, annotation, checkpoint hashing, and final checksums.

The template is tuned conservatively for a 140-thread/500-GB server: eight
library workers, 10 FastQC threads, eight Trim Galore threads, and eight
Bowtie2 plus four samtools threads produce a largest declared CPU request of
96. This leaves headroom for compression, filesystem work, and other users;
preflight reports the exact arithmetic for every run.

`RUN_FEATURE_ANNOTATION_SUMMARY=true` annotates every successful per-sample
caller peak set and optional primary consensus set. The output contains
exclusive feature counts/fractions/base-pair fractions and horizontal stacked
plots. Enhancer classification requires the matching `CCRE_BED_<GENOME>`;
without it, enhancer status is explicitly `not_evaluated` rather than silently
calling all distal peaks intergenic.

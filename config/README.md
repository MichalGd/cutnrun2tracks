# Configuration files

Copy `config.conf.template` to an untracked `config.conf`, then set every input,
output, and reference path used by the selected genome. The parser accepts only
the keys present in the template plus assembly-specific reference keys matching
`INDEX_*`, `FASTA_*`, `CHROM_SIZES_*`, `CANONICAL_CONTIGS_*`, `GTF_*`, `BLACKLIST_*`, `TSS_BED_*`, or
`CCRE_BED_*`.

`RUN_METAGENE=true` additionally requires `METAGENE_GENE_SET_MANIFEST` to name
a tab-separated gene-set manifest prepared with the utilities under
`common/metagene/`. The plotting stage consumes existing CPM or spike-in
bigWigs and does not perform normalization.

The examples contain illustrative `/data` and `/refs` paths and are not directly
executable. `cutrun_se.csv` requires `PEAK_CALLERS=macs3`.

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

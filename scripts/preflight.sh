#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

required=(python3 bowtie2 samtools bedtools "${PICARD_COMMAND}" "${MACS3_COMMAND}" bamCoverage Rscript)
needs_epic2=false
if [[ ",$PEAK_CALLERS," == *,epic2,* ]] && \
        awk -F '\t' 'NR>1 && $13=="FALSE" && ($9=="broad" || $9=="mixed") {found=1} END{exit !found}' "$SAMPLE_MANIFEST"; then
    needs_epic2=true
fi
is_true "$TRIM_ADAPTERS" && required+=(trim_galore)
is_true "$RUN_FASTQC" && required+=(fastqc)
is_true "$RUN_MULTIQC" && required+=(multiqc)
is_true "$GENERATE_COVERAGE_BIGWIGS" && required+=(bedGraphToBigWig)
[[ ",$PEAK_CALLERS," == *,seacr,* ]] && required+=("${SEACR_COMMAND}")
is_true "$needs_epic2" && required+=("${EPIC2_COMMAND}")
is_true "$RUN_PRESEQ" && required+=(preseq)
is_true "$RUN_CROSS_CORRELATION" && required+=("${PHANTOMPEAK_COMMAND}" gzip)
is_true "$RUN_REPLICATE_CORRELATION" && required+=(multiBamSummary plotCorrelation plotPCA)
is_true "$RUN_ATAQV_QC" && required+=(ataqv)
is_true "$GENERATE_ATAQV_VIEWER" && required+=(mkarv)
is_true "$RUN_TSS_SIGNAL_PROFILE" && required+=(computeMatrix plotProfile)
is_true "$RUN_METAGENE" && required+=(computeMatrix plotProfile plotHeatmap)
required+=(plotFingerprint)
is_true "$RUN_FEATURE_ANNOTATION_SUMMARY" && required+=(python3)
is_true "$RUN_MOTIF_ENRICHMENT" && die "RUN_MOTIF_ENRICHMENT is not implemented in this release"

missing=()
for command_name in "${required[@]}"; do
    if [[ "$command_name" == */* ]]; then
        [[ -x "$command_name" || -f "$command_name" ]] || missing+=("$command_name")
    elif ! command -v "$command_name" >/dev/null 2>&1; then
        missing+=("$command_name")
    fi
done
(( ${#missing[@]} == 0 )) || die "missing required commands: ${missing[*]}"

if is_true "$needs_epic2"; then
    epic2_probe_log="${OUTPUT_DIR}/00_metadata/epic2_preflight.log"
    "$EPIC2_COMMAND" --help >"$epic2_probe_log" 2>&1 || die "epic2 executable smoke check failed; see $epic2_probe_log"
    grep -q -- '--guess-bampe' "$epic2_probe_log" || die "epic2 executable does not expose --guess-bampe"
fi

visible_cpus="$(command -v nproc >/dev/null 2>&1 && nproc || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
[[ "$visible_cpus" =~ ^[1-9][0-9]*$ ]] || visible_cpus=1
[[ "$TOTAL_CPU_BUDGET" == "auto" ]] && cpu_budget="$visible_cpus" || cpu_budget="$TOTAL_CPU_BUDGET"
resource_table="${OUTPUT_DIR}/00_metadata/resource_budget.tsv"
printf 'stage\tparallel_jobs\tthreads_per_job\trequested_cpus\tbudget_cpus\tstatus\n' > "$resource_table"
resource_overcommit=0
check_stage_budget() {
    local stage="$1" jobs="$2" threads="$3" requested status=OK
    requested=$((jobs * threads))
    if (( requested > cpu_budget )); then
        status=OVERCOMMITTED; resource_overcommit=1
        warn "$stage requests up to $requested CPUs but TOTAL_CPU_BUDGET resolves to $cpu_budget"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$jobs" "$threads" "$requested" "$cpu_budget" "$status" >> "$resource_table"
}
(( THREADS_FASTQC > THREADS_TRIMGALORE )) && preprocess_threads="$THREADS_FASTQC" || preprocess_threads="$THREADS_TRIMGALORE"
alignment_threads=$((THREADS_BOWTIE2 + THREADS_SAMTOOLS))
(( THREADS_SAMTOOLS > THREADS_BAMCOVERAGE )) && qc_threads="$THREADS_SAMTOOLS" || qc_threads="$THREADS_BAMCOVERAGE"
check_stage_budget preprocess "$QC_SAMPLE_PARALLEL_JOBS" "$preprocess_threads"
check_stage_budget alignment "$THREADS_PARALLEL_JOBS" "$alignment_threads"
check_stage_budget filtering "$THREADS_PARALLEL_JOBS" "$THREADS_SAMTOOLS"
check_stage_budget cpm "$TRACK_PARALLEL_JOBS" "$THREADS_BAMCOVERAGE"
check_stage_budget peakcalling "$PEAKCALL_PARALLEL_JOBS" 1
check_stage_budget consensus "$MERGE_PARALLEL_JOBS" 1
check_stage_budget qc "$QC_SAMPLE_PARALLEL_JOBS" "$qc_threads"
check_stage_budget normalized_tracks "$NORMALIZED_TRACK_PARALLEL_JOBS" "$THREADS_BAMCOVERAGE"
check_stage_budget differential "$DIFFERENTIAL_PARALLEL_JOBS" 1
check_stage_budget annotation "$ANNOTATION_PARALLEL_JOBS" 1
[[ "$SPIKEIN_MODE" == "none" ]] || check_stage_budget spikein "$SPIKEIN_PARALLEL_JOBS" "$THREADS_BAMCOVERAGE"
is_true "$RUN_METAGENE" && check_stage_budget metagene "$METAGENE_PARALLEL_JOBS" "$METAGENE_THREADS_COMPUTEMATRIX"
if (( resource_overcommit )) && [[ "$RESOURCE_CHECK_MODE" == "fail" ]]; then
    die "configured parallelism exceeds TOTAL_CPU_BUDGET"
fi

if is_true "$RUN_METAGENE"; then
    [[ -s "$METAGENE_GENE_SET_MANIFEST" ]] || die "metagene gene-set manifest missing: $METAGENE_GENE_SET_MANIFEST"
    python3 -c 'import pyBigWig' || die "RUN_METAGENE requires the pyBigWig Python package"
fi
if is_true "$RUN_FEATURE_ANNOTATION_SUMMARY"; then
    python3 -c 'import matplotlib' || die "RUN_FEATURE_ANNOTATION_SUMMARY requires matplotlib"
fi

mkdir -p "${OUTPUT_DIR}/00_metadata"
versions="${OUTPUT_DIR}/00_metadata/software_versions.tsv"
printf 'tool\tversion\n' > "$versions"
version_commands=(python3 bowtie2 samtools bedtools bamCoverage Rscript "$MACS3_COMMAND")
is_true "$needs_epic2" && version_commands+=("$EPIC2_COMMAND")
is_true "$RUN_CROSS_CORRELATION" && version_commands+=("$PHANTOMPEAK_COMMAND")
is_true "$RUN_METAGENE" && version_commands+=(computeMatrix plotProfile plotHeatmap)
for command_name in "${version_commands[@]}"; do
    version="$($command_name --version 2>&1 | head -n 1 || true)"
    printf '%s\t%s\n' "$command_name" "${version//$'\t'/ }" >> "$versions"
done

tail -n +2 "$SAMPLE_MANIFEST" | while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    index="$(reference_value INDEX "$genome")"
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    canonical_contigs="$(reference_value CANONICAL_CONTIGS "$genome")"
    fasta="$(reference_value FASTA "$genome")"
    gtf="$(reference_value GTF "$genome")"
    [[ -s "$chrom_sizes" && -s "$canonical_contigs" && -s "$fasta" && -s "$gtf" && -s "$blacklist" ]] || \
        die "missing reference file for $sample_key"
    if is_true "$RUN_CCRE_ANNOTATION"; then
        ccre="$(optional_reference_value CCRE_BED "$genome")"
        [[ -s "$ccre" ]] || die "RUN_CCRE_ANNOTATION requires CCRE_BED for $genome"
    fi
    [[ -f "${index}.1.bt2" || -f "${index}.1.bt2l" ]] || die "Bowtie2 index prefix invalid: $index"
    if is_true "$needs_epic2"; then
        effective="$(reference_value EFFECTIVE_GENOME_SIZE "$genome")"
        [[ "$effective" =~ ^[1-9][0-9]*$ ]] || die "effective genome size must be a positive integer for $genome"
    fi
done

if [[ "$SPIKEIN_MODE" != "none" ]]; then
    [[ -s "$SPIKEIN_FASTA" && -s "$SPIKEIN_CHROM_SIZES" && -s "$SPIKEIN_ALLOWED_CONTIGS" ]] || \
        die "spike-in reference manifest files are incomplete"
    [[ -f "${SPIKEIN_INDEX}.1.bt2" || -f "${SPIKEIN_INDEX}.1.bt2l" ]] || \
        die "composite Bowtie2 index prefix invalid: $SPIKEIN_INDEX"
fi

profiles="$(awk -F '\t' 'NR>1 {seen[$6]=1} END {for (value in seen) print value}' "$SAMPLE_MANIFEST" | LC_ALL=C sort | paste -sd, -)"
printf 'status\tprofiles\tspikein_mode\tpeak_callers\tcpu_budget\nSUCCESS\t%s\t%s\t%s\t%s\n' \
    "$profiles" "$SPIKEIN_MODE" "$PEAK_CALLERS" "$cpu_budget" > "${OUTPUT_DIR}/00_metadata/preflight_status.tsv"
note "Preflight successful"

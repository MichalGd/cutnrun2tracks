#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/01_fastq_qc/raw" "${OUTPUT_DIR}/01_fastq_qc/trimmed" \
    "${OUTPUT_DIR}/01_fastq_qc/multiqc" "${OUTPUT_DIR}/02_trimmed_fastq" "${OUTPUT_DIR}/logs/preprocess"

merge_fastqs() {
    local list="${1:?list}" output="${2:?output}" path
    : > "$output"
    IFS=';' read -r -a paths <<< "$list"
    for path in "${paths[@]}"; do
        [[ -s "$path" ]] || die "FASTQ missing: $path"
        gzip -t "$path"
        # Concatenated gzip members are valid FASTQ.gz and avoid recompression.
        dd if="$path" of="$output" oflag=append conv=notrunc status=none
    done
    gzip -t "$output"
}

worker() {
    local key="$1" layout="$2" fq1_list="$3" fq2_list="$4"
    local work="${OUTPUT_DIR}/02_trimmed_fastq"
    local raw1="${work}/${key}.R1.merged.fastq.gz"
    local raw2="${work}/${key}.R2.merged.fastq.gz" log="${OUTPUT_DIR}/logs/preprocess/${key}.log"
    merge_fastqs "$fq1_list" "$raw1"
    [[ "$layout" == "PE" ]] && merge_fastqs "$fq2_list" "$raw2"
    if is_true "$RUN_FASTQC"; then
        if [[ "$layout" == "PE" ]]; then
            run_logged fastqc --threads "$THREADS_FASTQC" --outdir "${OUTPUT_DIR}/01_fastq_qc/raw" "$raw1" "$raw2" >"$log" 2>&1
        else
            run_logged fastqc --threads "$THREADS_FASTQC" --outdir "${OUTPUT_DIR}/01_fastq_qc/raw" "$raw1" >"$log" 2>&1
        fi
    fi
    if is_true "$TRIM_ADAPTERS"; then
        local extra=()
        [[ -n "$TRIMGALORE_EXTRA_ARGS" ]] && read -r -a extra <<< "$TRIMGALORE_EXTRA_ARGS"
        if [[ "$layout" == "PE" ]]; then
            run_logged trim_galore --paired --cores "$THREADS_TRIMGALORE" --length "$MIN_TRIMMED_LENGTH" --output_dir "$work" \
                "${extra[@]}" "$raw1" "$raw2" >>"$log" 2>&1
            mv "${work}/${key}.R1.merged_val_1.fq.gz" "${work}/${key}.R1.trimmed.fastq.gz"
            mv "${work}/${key}.R2.merged_val_2.fq.gz" "${work}/${key}.R2.trimmed.fastq.gz"
        else
            run_logged trim_galore --cores "$THREADS_TRIMGALORE" --length "$MIN_TRIMMED_LENGTH" --output_dir "$work" \
                "${extra[@]}" "$raw1" >>"$log" 2>&1
            mv "${work}/${key}.R1.merged_trimmed.fq.gz" "${work}/${key}.R1.trimmed.fastq.gz"
        fi
    else
        ln -sfn "$(basename "$raw1")" "${work}/${key}.R1.trimmed.fastq.gz"
        [[ "$layout" == "PE" ]] && ln -sfn "$(basename "$raw2")" "${work}/${key}.R2.trimmed.fastq.gz"
    fi
    if is_true "$RUN_FASTQC"; then
        local trimmed1="${work}/${key}.R1.trimmed.fastq.gz" trimmed2=""
        [[ "$layout" == "PE" ]] && trimmed2="${work}/${key}.R2.trimmed.fastq.gz"
        if [[ "$layout" == "PE" ]]; then
            run_logged fastqc --threads "$THREADS_FASTQC" --outdir "${OUTPUT_DIR}/01_fastq_qc/trimmed" "$trimmed1" "$trimmed2" >>"$log" 2>&1
        else
            run_logged fastqc --threads "$THREADS_FASTQC" --outdir "${OUTPUT_DIR}/01_fastq_qc/trimmed" "$trimmed1" >>"$log" 2>&1
        fi
    fi
}

parallel_require_positive_integer "THREADS_FASTQC" "$THREADS_FASTQC"
parallel_require_positive_integer "THREADS_TRIMGALORE" "$THREADS_TRIMGALORE"
parallel_pool_init "$QC_SAMPLE_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" ]] && continue
    parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$fastq_1_list" "$fastq_2_list"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all

if is_true "$RUN_MULTIQC"; then
    run_logged multiqc --force --outdir "${OUTPUT_DIR}/01_fastq_qc/multiqc" "${OUTPUT_DIR}/01_fastq_qc" \
        >"${OUTPUT_DIR}/logs/preprocess/multiqc.log" 2>&1
fi
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/02_trimmed_fastq/stage_status.tsv"

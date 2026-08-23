#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/05_peaks/per_sample" "${OUTPUT_DIR}/logs/peakcalling"

genome_size() {
    local genome="$1" chrom_sizes
    if [[ "$MACS3_GENOME_SIZE" != "auto" ]]; then
        printf '%s' "$MACS3_GENOME_SIZE"
        return
    fi
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    awk '{total += $2} END {printf "%.0f", total}' "$chrom_sizes"
}

fragment_bedgraph() {
    local bam="$1" output="$2" tmp
    tmp="$(mktemp -d "$(dirname "$output")/.seacr-fragments.XXXXXX")"
    samtools view -h -f 2 -F 2828 "$bam" |
        awk -v maximum="$SEACR_MAX_FRAGMENT" 'BEGIN{OFS="\t"} /^@/ {print; next} {t=$9; if(t<0)t=-t; if(t>0 && t<=maximum)print}' |
        samtools view -b -o "$tmp/fragments.bam" -
    samtools sort -@ "$THREADS_SAMTOOLS" -o "$tmp/sorted.bam" "$tmp/fragments.bam"
    bedtools genomecov -bg -pc -ibam "$tmp/sorted.bam" | awk '$4>0' > "$output"
    [[ -s "$output" ]] || die "SEACR fragment bedGraph is empty: $bam"
    rm -rf "$tmp"
}

call_macs3() {
    local key="$1" layout="$2" genome="$3" target_class="$4" target_bam="$5" control_bam="$6" out="$7"
    mkdir -p "$out"
    local format=() control=() common=()
    [[ "$control_bam" != "." ]] && control=(-c "$control_bam")
    if [[ "$layout" == "PE" ]]; then
        format=(-f BAMPE)
    elif [[ "$ASSAY_PROFILE" == "cutrun" ]]; then
        format=(-f BAM --nomodel --shift "$MACS3_CUTRUN_SE_SHIFT" --extsize "$MACS3_CUTRUN_SE_EXTSIZE")
    else
        format=(-f BAM --nomodel --shift "$MACS3_CUTTAG_SE_SHIFT" --extsize "$MACS3_CUTTAG_SE_EXTSIZE")
    fi
    common=(callpeak -t "$target_bam" "${control[@]}" "${format[@]}" -g "$(genome_size "$genome")"
        -q "$MACS3_QVALUE" --keep-dup "$MACS3_KEEP_DUP" --outdir "$out")
    if [[ "$target_class" == "narrow" || "$target_class" == "mixed" ]]; then
        local summit_args=()
        is_true "$MACS3_CALL_SUMMITS" && summit_args+=(--call-summits)
        run_logged "$MACS3_COMMAND" "${common[@]}" -n "${key}.macs3.narrow" "${summit_args[@]}" \
            >"${OUTPUT_DIR}/logs/peakcalling/${key}.macs3.narrow.log" 2>&1
        local narrow="${out}/${key}.macs3.narrow_peaks.narrowPeak"
        if [[ -s "$narrow" ]]; then cut -f1-3 "$narrow" > "${out}/${key}.macs3.narrow.bed"; fi
    fi
    if [[ "$target_class" == "broad" || "$target_class" == "mixed" ]]; then
        run_logged "$MACS3_COMMAND" "${common[@]}" -n "${key}.macs3.broad" --broad --broad-cutoff "$MACS3_BROAD_CUTOFF" \
            >"${OUTPUT_DIR}/logs/peakcalling/${key}.macs3.broad.log" 2>&1
        local broad="${out}/${key}.macs3.broad_peaks.broadPeak"
        if [[ -s "$broad" ]]; then cut -f1-3 "$broad" > "${out}/${key}.macs3.broad.bed"; fi
    fi
}

call_seacr() {
    local key="$1" target_bam="$2" control_bam="$3" out="$4"
    mkdir -p "$out"
    local target_bg="${out}/${key}.target.fragments.bedGraph" control_arg prefix produced
    fragment_bedgraph "$target_bam" "$target_bg"
    if [[ "$control_bam" != "." ]]; then
        control_arg="${out}/${key}.control.fragments.bedGraph"
        fragment_bedgraph "$control_bam" "$control_arg"
    else
        control_arg="$SEACR_NO_CONTROL_THRESHOLD"
    fi
    prefix="${out}/${key}.seacr"
    if [[ -f "$SEACR_COMMAND" ]]; then
        run_logged bash "$SEACR_COMMAND" "$target_bg" "$control_arg" "$SEACR_CONTROL_NORMALIZATION" "$SEACR_MODE" "$prefix" \
            >"${OUTPUT_DIR}/logs/peakcalling/${key}.seacr.log" 2>&1
    else
        run_logged "$SEACR_COMMAND" "$target_bg" "$control_arg" "$SEACR_CONTROL_NORMALIZATION" "$SEACR_MODE" "$prefix" \
            >"${OUTPUT_DIR}/logs/peakcalling/${key}.seacr.log" 2>&1
    fi
    produced="${prefix}.${SEACR_MODE}.bed"
    if [[ -s "$produced" ]]; then
        cut -f1-3 "$produced" > "${out}/${key}.seacr.narrow.bed"
    elif ! is_true "$ALLOW_EMPTY_PEAKS"; then
        die "SEACR produced no peaks for $key"
    else
        : > "${out}/${key}.seacr.narrow.bed"
    fi
}

worker() {
    local key="$1" layout="$2" genome="$3" target_class="$4" control_key="$5" primary_caller="$6" primary_class="$7"
    local target_bam control_bam="." root
    target_bam="$(analysis_bam_path "$key")"
    [[ "$control_key" != "." ]] && control_bam="$(analysis_bam_path "$control_key")"
    root="${OUTPUT_DIR}/05_peaks/per_sample/${key}"
    if [[ ",$PEAK_CALLERS," == *,macs3,* ]]; then
        call_macs3 "$key" "$layout" "$genome" "$target_class" "$target_bam" "$control_bam" "${root}/macs3"
    fi
    if [[ ",$PEAK_CALLERS," == *,seacr,* ]]; then
        [[ "$layout" == "PE" ]] || die "SEACR requested for SE sample $key"
        call_seacr "$key" "$target_bam" "$control_bam" "${root}/seacr"
    fi
    local primary="${root}/${primary_caller}/${key}.${primary_caller}.${primary_class}.bed"
    if [[ ! -s "$primary" ]]; then
        if is_true "$ALLOW_EMPTY_PEAKS"; then
            : > "$primary"
        else
            die "primary peak file missing or empty for $key: $primary"
        fi
    fi
    printf 'sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\n%s\t%s\t%s\t%s\tSUCCESS\n' \
        "$key" "$control_key" "$primary_caller" "$primary_class" > "${root}/peakcall_metadata.tsv"
}

parallel_pool_init "$PEAKCALL_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" || "$is_control" == "TRUE" ]] && continue
    parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome" "$target_class" \
        "$control_key" "$primary_caller" "$primary_class"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/05_peaks/per_sample/stage_status.tsv"

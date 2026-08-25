#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

NORMALIZATION_LOG_DIR="${OUTPUT_DIR}/logs/normalized_tracks"
FAMILY_STATUS="${OUTPUT_DIR}/04_tracks/normalized_track_family_status.tsv"
normalization_skips=0

mkdir -p "$NORMALIZATION_LOG_DIR" "${OUTPUT_DIR}/04_tracks/deseq2_consensus" \
    "${OUTPUT_DIR}/04_tracks/deseq2_robust_cpm"
rm -f -- "${OUTPUT_DIR}/04_tracks/stage_status.tsv"
printf 'cohort_id\tpolicy\tstatus\treason\tlog\n' > "$FAMILY_STATUS"

invalidate_factor_outputs() {
    local tables="$1"
    rm -f -- \
        "${tables}/normalization_factors.tsv" \
        "${tables}/raw_counts.tsv.gz" \
        "${tables}/normalized_counts.tsv.gz" \
        "${tables}/session_info.txt" \
        "${tables}/consensus_count_sums.tsv"
}

record_family_status() {
    local cohort="$1" policy="$2" status="$3" reason="$4" log_path="$5"
    printf '%s\t%s\t%s\t%s\t%s\n' "$cohort" "$policy" "$status" "$reason" "$log_path" >> "$FAMILY_STATUS"
}

skip_or_fail_family() {
    local cohort="$1" policy="$2" family_dir="$3" reason="$4" log_path="$5"
    normalization_skips=$((normalization_skips + 1))
    if is_true "$REQUIRE_ALL_ENABLED_TRACKS"; then
        printf '{"status":"FAILED","reason":"%s","cohort_id":"%s","policy":"%s","log":"%s"}\n' \
            "$reason" "$cohort" "$policy" "$log_path" > "${family_dir}/FAILED.json"
        record_family_status "$cohort" "$policy" FAILED "$reason" "$log_path"
        printf 'status\tskipped_families\nFAILED\t%s\n' "$normalization_skips" \
            > "${OUTPUT_DIR}/04_tracks/stage_status.tsv"
        return 1
    fi
    printf '{"status":"SKIPPED","reason":"%s","cohort_id":"%s","policy":"%s","log":"%s"}\n' \
        "$reason" "$cohort" "$policy" "$log_path" > "${family_dir}/SKIPPED.json"
    record_family_status "$cohort" "$policy" SKIPPED "$reason" "$log_path"
    warn "skipping normalized-track family for cohort $cohort ($policy): $reason${log_path:+; see $log_path}"
    return 0
}

track_from_factor() {
    local sample_key="$1" layout="$2" genome="$3" bam="$4" scale="$5" output_stem="$6"
    local chrom_sizes tmp sorted extra=() args=()
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    mkdir -p "$(dirname "$output_stem")"
    tmp="${output_stem}.unsorted.bedGraph"
    sorted="${output_stem}.bedGraph"
    args=(--bam "$bam" --outFileName "$tmp" --outFileFormat bedgraph --binSize "$TRACK_BIN_SIZE"
        --scaleFactor "$scale" --numberOfProcessors "$THREADS_BAMCOVERAGE")
    [[ "$layout" == "PE" ]] && args+=(--samFlagInclude 66 --extendReads)
    [[ -n "$BAMCOVERAGE_COMMON_ARGS" ]] && read -r -a extra <<< "$BAMCOVERAGE_COMMON_ARGS"
    run_logged bamCoverage "${args[@]}" "${extra[@]}"
    bedtools sort -faidx "$chrom_sizes" -i "$tmp" > "$sorted"
    awk 'BEGIN{ok=1} NF!=4 || $2<0 || $3<=$2 || $4<0 || $4!=$4 {ok=0} END{exit !ok}' "$sorted" || \
        die "invalid normalized bedGraph: $sorted"
    is_true "$GENERATE_COVERAGE_BIGWIGS" && run_logged bedGraphToBigWig "$sorted" "$chrom_sizes" "${output_stem}.bw"
    if is_true "$GENERATE_COVERAGE_BEDGRAPHS"; then gzip -c "$sorted" > "${output_stem}.bedGraph.gz"; fi
    is_true "$KEEP_RAW_BEDGRAPH" || rm -f "$tmp" "$sorted"
}

run_family() {
    local cohort="$1" genome="$2" layout="$3" caller="$4" peak_class="$5" policy="$6" family_dir="$7" label="$8" emit_tracks="$9"
    local consensus_dir consensus tables factors factor_log
    consensus_dir="${OUTPUT_DIR}/05_peaks/consensus/${cohort}/${caller}/${peak_class}"
    consensus="$(find "$consensus_dir" -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
    tables="${family_dir}/tables"
    factor_log="${NORMALIZATION_LOG_DIR}/${cohort}.${policy}.factors.log"
    mkdir -p "$family_dir" "$tables"
    rm -f -- "${family_dir}/SKIPPED.json" "${family_dir}/FAILED.json"
    invalidate_factor_outputs "$tables"
    if [[ -z "$consensus" ]]; then
        skip_or_fail_family "$cohort" "$policy" "$family_dir" "consensus unavailable" "."
        return $?
    fi
    if ! run_logged Rscript "${SCRIPT_DIR}/consensus_track_factors.R" "$SAMPLE_MANIFEST" "$cohort" "$consensus" \
        "$OUTPUT_DIR" "$tables" "$policy" > "$factor_log" 2>&1; then
        skip_or_fail_family "$cohort" "$policy" "$family_dir" "normalization factor calculation failed" "$factor_log"
        return $?
    fi
    factors="${tables}/normalization_factors.tsv"
    if [[ ! -s "$factors" ]]; then
        skip_or_fail_family "$cohort" "$policy" "$family_dir" "normalization factor table is missing or empty" "$factor_log"
        return $?
    fi
    while IFS=$'\t' read -r sample_key factor_cohort factor_policy signal_unit size_factor consensus_scale count_sum geometric_mean effective robust_scale; do
        [[ "$sample_key" == "sample_key" ]] && continue
        local scale bam stem
        if [[ "$policy" == "analysis" ]]; then scale="$consensus_scale"; else scale="$robust_scale"; fi
        bam="$(policy_bam_path "$sample_key" "$policy")"
        stem="${family_dir}/${sample_key}.${label}"
        is_true "$emit_tracks" && track_from_factor "$sample_key" "$layout" "$genome" "$bam" "$scale" "$stem"
    done < "$factors"
    record_family_status "$cohort" "$policy" SUCCESS . "$factor_log"
}

while IFS=$'\t' read -r cohort cohort_key genome assay_profile factor antibody_id layout target_class duplicate_policy caller peak_class n_samples sample_keys conditions; do
    [[ "$cohort" == "cohort_id" ]] && continue
    if is_true "$GENERATE_DESEQ2_CONSENSUS_TRACKS" || is_true "$RUN_DIFFBIND" || is_true "$RUN_DESEQ2_ENRICHMENT"; then
        run_family "$cohort" "$genome" "$layout" "$caller" "$peak_class" analysis \
            "${OUTPUT_DIR}/04_tracks/deseq2_consensus/${cohort}" DESeq2Consensus "$GENERATE_DESEQ2_CONSENSUS_TRACKS"
    fi
    if is_true "$GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS"; then
        run_family "$cohort" "$genome" "$layout" "$caller" "$peak_class" permissive \
            "${OUTPUT_DIR}/04_tracks/deseq2_robust_cpm/permissive/${cohort}" DESeq2RobustCPM.Permissive true
    fi
    if is_true "$GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS"; then
        run_family "$cohort" "$genome" "$layout" "$caller" "$peak_class" intermediate \
            "${OUTPUT_DIR}/04_tracks/deseq2_robust_cpm/intermediate/${cohort}" DESeq2RobustCPM.Intermediate true
    fi
    if is_true "$GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS"; then
        run_family "$cohort" "$genome" "$layout" "$caller" "$peak_class" stringent \
            "${OUTPUT_DIR}/04_tracks/deseq2_robust_cpm/stringent/${cohort}" DESeq2RobustCPM.Stringent true
    fi
done < "$COHORT_MANIFEST"
if (( normalization_skips > 0 )); then
    printf 'status\tskipped_families\nCOMPLETED_WITH_WARNINGS\t%s\n' "$normalization_skips" \
        > "${OUTPUT_DIR}/04_tracks/stage_status.tsv"
else
    printf 'status\tskipped_families\nSUCCESS\t0\n' > "${OUTPUT_DIR}/04_tracks/stage_status.tsv"
fi

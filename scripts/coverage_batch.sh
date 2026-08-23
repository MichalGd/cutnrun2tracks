#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/04_tracks/cpm" "${OUTPUT_DIR}/logs/coverage"

worker() {
    local key="$1" layout="$2" genome="$3"
    local bam count scale chrom_sizes tmp sorted bedgraph bigwig
    bam="$(analysis_bam_path "$key")"
    count="$(signal_count "$bam" "$layout")"
    (( count > 0 )) || die "zero analysis observations for $key"
    scale="$(awk -v n="$count" 'BEGIN {printf "%.15g", 1000000/n}')"
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    tmp="${OUTPUT_DIR}/04_tracks/cpm/${key}.CPM.unsorted.bedGraph"
    sorted="${OUTPUT_DIR}/04_tracks/cpm/${key}.CPM.bedGraph"
    bedgraph="${sorted}.gz"
    bigwig="${OUTPUT_DIR}/04_tracks/cpm/${key}.CPM.bw"
    local args=(--bam "$bam" --outFileName "$tmp" --outFileFormat bedgraph --binSize "$TRACK_BIN_SIZE"
        --scaleFactor "$scale" --numberOfProcessors "$THREADS_BAMCOVERAGE")
    [[ "$layout" == "PE" ]] && args+=(--samFlagInclude 66 --extendReads)
    local extra=()
    [[ -n "$BAMCOVERAGE_COMMON_ARGS" ]] && read -r -a extra <<< "$BAMCOVERAGE_COMMON_ARGS"
    run_logged bamCoverage "${args[@]}" "${extra[@]}" >"${OUTPUT_DIR}/logs/coverage/${key}.log" 2>&1
    bedtools sort -faidx "$chrom_sizes" -i "$tmp" > "$sorted"
    awk 'BEGIN{ok=1} NF!=4 || $2<0 || $3<=$2 || $4<0 || $4!=$4 {ok=0} END{exit !ok}' "$sorted" || \
        die "invalid bedGraph: $sorted"
    if is_true "$GENERATE_COVERAGE_BIGWIGS"; then
        run_logged bedGraphToBigWig "$sorted" "$chrom_sizes" "$bigwig"
    fi
    if is_true "$GENERATE_COVERAGE_BEDGRAPHS"; then
        gzip -c "$sorted" > "$bedgraph"
        gzip -t "$bedgraph"
    fi
    is_true "$KEEP_RAW_BEDGRAPH" || rm -f "$tmp" "$sorted"
    printf 'sample_key\tsignal_unit\tsignal_count\tscale\tformula\n%s\t%s\t%s\t%s\tC*1e6/L\n' \
        "$key" "$([[ "$layout" == "PE" ]] && echo fragment || echo read)" "$count" "$scale" \
        > "${OUTPUT_DIR}/04_tracks/cpm/${key}.normalization_metadata.tsv"
}

if is_true "$GENERATE_CPM_TRACKS"; then
    parallel_pool_init "$TRACK_PARALLEL_JOBS"
    while IFS=$'\t' read -r sample_key sample_id replicate layout genome rest; do
        [[ "$sample_key" == "sample_key" ]] && continue
        parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome"
    done < "$SAMPLE_MANIFEST"
    parallel_pool_wait_all
else
    printf '{"status":"SKIPPED","reason":"GENERATE_CPM_TRACKS=false"}\n' > "${OUTPUT_DIR}/04_tracks/cpm/SKIPPED.json"
fi
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/04_tracks/cpm/stage_status.tsv"

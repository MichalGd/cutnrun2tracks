#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

for branch in q0_dup-retained q0_dup-removed q30_dup-retained q30_dup-removed; do
    mkdir -p "${OUTPUT_DIR}/03_alignment/filtered/${branch}"
done
mkdir -p "${OUTPUT_DIR}/03_alignment/marked" "${OUTPUT_DIR}/03_alignment/analysis" \
    "${OUTPUT_DIR}/03_alignment/metrics" "${OUTPUT_DIR}/logs/filtering"

filter_branch() (
    local marked="$1" output="$2" layout="$3" genome="$4" blacklist="$5" mapq="$6" remove_duplicates="$7"
    local tmp canonical_file include=0 exclude contigs=()
    tmp="$(mktemp -d "${OUTPUT_DIR}/03_alignment/filtered/.filter.XXXXXX")"
    trap 'rm -rf -- "$tmp"' EXIT
    rm -f "$output" "${output}.bai" "${output%.bam}.bai"
    if [[ "$layout" == "PE" ]]; then
        include=2
        [[ "$remove_duplicates" == "true" ]] && exclude=3852 || exclude=2828
    else
        [[ "$remove_duplicates" == "true" ]] && exclude=3844 || exclude=2820
    fi
    canonical_file="$(reference_value CANONICAL_CONTIGS "$genome")"
    mapfile -t contigs < <(awk -v remove_mito="$REMOVE_MITO" '
        NF>=2 && !(remove_mito=="true" && ($1=="chrM" || $1=="MT" || $1=="M")) {print $1}
        NF==1 && !(remove_mito=="true" && ($1=="chrM" || $1=="MT" || $1=="M")) {print $1}
    ' "$canonical_file")
    (( ${#contigs[@]} > 0 )) || die "canonical contig list is empty: $canonical_file"
    # Region selection requires an indexed input. Apply flags, MAPQ, and
    # canonical-contig selection in one pass from the already indexed marked BAM.
    samtools view -@ "$THREADS_SAMTOOLS" -b -q "$mapq" -f "$include" -F "$exclude" \
        -o "$tmp/canonical.bam" "$marked" "${contigs[@]}"
    bedtools intersect -v -abam "$tmp/canonical.bam" -b "$blacklist" > "$tmp/no_blacklist.bam"
    if [[ "$layout" == "PE" ]]; then
        samtools sort -n -@ "$THREADS_SAMTOOLS" -o "$tmp/name.bam" "$tmp/no_blacklist.bam"
        samtools fixmate -r -@ "$THREADS_SAMTOOLS" "$tmp/name.bam" "$tmp/fixmate.bam"
        samtools view -@ "$THREADS_SAMTOOLS" -b -f 2 -F "$exclude" -o "$tmp/paired.bam" "$tmp/fixmate.bam"
        samtools sort -@ "$THREADS_SAMTOOLS" -o "$output" "$tmp/paired.bam"
    else
        samtools sort -@ "$THREADS_SAMTOOLS" -o "$output" "$tmp/no_blacklist.bam"
    fi
    samtools quickcheck "$output"
    samtools index -@ "$THREADS_SAMTOOLS" "$output"
    local count
    count="$(signal_count "$output" "$layout")"
    if (( count == 0 )) && ! is_true "$ALLOW_EMPTY_FILTERED_BAM"; then
        die "filtering removed every signal unit: $output"
    fi
)

worker() {
    local key="$1" layout="$2" genome="$3" duplicate_policy="$4" blacklist="$5"
    local sorted="${OUTPUT_DIR}/03_alignment/sorted/${key}.host.sorted.bam"
    local marked="${OUTPUT_DIR}/03_alignment/marked/${key}.host.marked.bam"
    local metrics="${OUTPUT_DIR}/03_alignment/metrics/${key}.duplicate_metrics.txt"
    [[ -s "$sorted" ]] || die "sorted BAM missing: $sorted"
    local picard_log="${OUTPUT_DIR}/logs/filtering/${key}.picard.log"
    if [[ -s "$marked" && -s "$metrics" ]] && \
        samtools quickcheck "$marked" && samtools idxstats "$marked" >/dev/null 2>&1; then
        note "Reusing validated marked BAM for $key" >> "$picard_log"
    else
        rm -f "$marked" "${marked}.bai" "${marked%.bam}.bai" "$metrics"
        run_logged "$PICARD_COMMAND" MarkDuplicates I="$sorted" O="$marked" M="$metrics" \
            REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
            >"$picard_log" 2>&1
    fi
    samtools quickcheck "$marked"
    samtools idxstats "$marked" >/dev/null
    local q0r="${OUTPUT_DIR}/03_alignment/filtered/q0_dup-retained/${key}.host.q0.dup-retained.bam"
    local q0d="${OUTPUT_DIR}/03_alignment/filtered/q0_dup-removed/${key}.host.q0.dup-removed.bam"
    local q30r="${OUTPUT_DIR}/03_alignment/filtered/q30_dup-retained/${key}.host.q30.dup-retained.bam"
    local q30d="${OUTPUT_DIR}/03_alignment/filtered/q30_dup-removed/${key}.host.q30.dup-removed.bam"
    filter_branch "$marked" "$q0r" "$layout" "$genome" "$blacklist" "$PERMISSIVE_MIN_MAPQ" false
    filter_branch "$marked" "$q0d" "$layout" "$genome" "$blacklist" "$INTERMEDIATE_MIN_MAPQ" true
    filter_branch "$marked" "$q30r" "$layout" "$genome" "$blacklist" "$MIN_MAPQ" false
    filter_branch "$marked" "$q30d" "$layout" "$genome" "$blacklist" "$MIN_MAPQ" true
    local analysis_source
    [[ "$duplicate_policy" == "retain" ]] && analysis_source="$q30r" || analysis_source="$q30d"
    ln -sfn "$analysis_source" "${OUTPUT_DIR}/03_alignment/analysis/${key}.host.analysis.bam"
    ln -sfn "${analysis_source}.bai" "${OUTPUT_DIR}/03_alignment/analysis/${key}.host.analysis.bam.bai"
    printf 'sample_key\tq0_dup_retained\tq0_dup_removed\tq30_dup_retained\tq30_dup_removed\tanalysis_policy\n' \
        > "${OUTPUT_DIR}/03_alignment/metrics/${key}.filter_counts.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$key" \
        "$(signal_count "$q0r" "$layout")" "$(signal_count "$q0d" "$layout")" \
        "$(signal_count "$q30r" "$layout")" "$(signal_count "$q30d" "$layout")" "$duplicate_policy" \
        >> "${OUTPUT_DIR}/03_alignment/metrics/${key}.filter_counts.tsv"
}

parallel_pool_init "$THREADS_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist rest; do
    [[ "$sample_key" == "sample_key" ]] && continue
    parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome" "$duplicate_policy" "$blacklist"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/03_alignment/analysis/stage_status.tsv"

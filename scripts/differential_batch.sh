#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

status_dir="${OUTPUT_DIR}/08_differential/.cohort_status"
mkdir -p "${OUTPUT_DIR}/08_differential" "${OUTPUT_DIR}/logs/differential" "$status_dir"
rm -f -- "$status_dir"/*.tsv 2>/dev/null || true

differential_worker() {
    local cohort="$1" caller="$2" peak_class="$3" failures=0 skips=0
    local consensus_dir consensus root counts upstream_skip spike_file
    consensus_dir="${OUTPUT_DIR}/05_peaks/consensus/${cohort}/${caller}/${peak_class}"
    consensus="$(find "$consensus_dir" -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
    root="${OUTPUT_DIR}/08_differential/${cohort}/${peak_class}"
    mkdir -p "$root/primary_target_only" "$root/sensitivity_control_subtracted" \
        "$root/sensitivity_target_control_interaction" "$root/concordance"
    rm -f -- "$root/SKIPPED.json" "$root/FAILED.json"
    if [[ -z "$consensus" ]]; then
        printf '{"status":"SKIPPED","reason":"consensus unavailable"}\n' > "$root/SKIPPED.json"
        printf '%s\t0\t1\n' "$cohort" > "$status_dir/${cohort}.tsv"
        return 0
    fi
    counts="${OUTPUT_DIR}/04_tracks/deseq2_consensus/${cohort}/tables/raw_counts.tsv.gz"
    if [[ ! -s "$counts" ]]; then
        upstream_skip="${OUTPUT_DIR}/04_tracks/deseq2_consensus/${cohort}/SKIPPED.json"
        if [[ -s "$upstream_skip" ]]; then
            printf '{"status":"SKIPPED","reason":"consensus normalization unavailable","upstream_status":"%s"}\n' \
                "$upstream_skip" > "$root/SKIPPED.json"
            skips=1
        else
            printf '{"status":"FAILED","reason":"raw consensus counts unexpectedly unavailable"}\n' > "$root/FAILED.json"
            failures=1
        fi
        printf '%s\t%s\t%s\n' "$cohort" "$failures" "$skips" > "$status_dir/${cohort}.tsv"
        return 0
    fi
    spike_file="."
    [[ "$DIFFERENTIAL_NORMALIZATION" == "spikein" ]] && spike_file="${OUTPUT_DIR}/06_qc/spikein/spikein_scaling.tsv"
    if is_true "$RUN_DESEQ2_ENRICHMENT"; then
        if ! run_logged Rscript "${SCRIPT_DIR}/deseq2_enrichment_analysis.R" "$SAMPLE_MANIFEST" "$cohort" "$counts" \
            "$root/primary_target_only/deseq2_enrichment" "$DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION" \
            "$DIFFERENTIAL_ALPHA" "$DIFFERENTIAL_MIN_ABS_LOG2FC" "${DIFFERENTIAL_BLOCK_COLUMNS:-.}" \
            "${DIFFERENTIAL_CONDITION_ORDER:-.}" "${DIFFERENTIAL_REFERENCE_CONDITION:-.}" "$spike_file" \
            >"${OUTPUT_DIR}/logs/differential/${cohort}.deseq2.log" 2>&1; then
            failures=$((failures+1)); mkdir -p "$root/primary_target_only/deseq2_enrichment"
            printf '{"status":"FAILED","module":"DESeq2Enrichment"}\n' > "$root/primary_target_only/deseq2_enrichment/FAILED.json"
        fi
    fi
    if is_true "$RUN_DIFFBIND"; then
        if ! run_logged Rscript "${SCRIPT_DIR}/diffbind_analysis.R" "$SAMPLE_MANIFEST" "$cohort" "$consensus" "$OUTPUT_DIR" \
            "$root/primary_target_only/diffbind" "$DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION" "$DIFFERENTIAL_ALPHA" \
            "${DIFFERENTIAL_BLOCK_COLUMNS:-.}" false >"${OUTPUT_DIR}/logs/differential/${cohort}.diffbind.log" 2>&1; then
            failures=$((failures+1)); mkdir -p "$root/primary_target_only/diffbind"
            printf '{"status":"FAILED","module":"DiffBind"}\n' > "$root/primary_target_only/diffbind/FAILED.json"
        fi
        if is_true "$RUN_CONTROL_SUBTRACTED_SENSITIVITY"; then
            if ! run_logged Rscript "${SCRIPT_DIR}/diffbind_analysis.R" "$SAMPLE_MANIFEST" "$cohort" "$consensus" "$OUTPUT_DIR" \
                "$root/sensitivity_control_subtracted/diffbind" "$DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION" "$DIFFERENTIAL_ALPHA" \
                "${DIFFERENTIAL_BLOCK_COLUMNS:-.}" true >"${OUTPUT_DIR}/logs/differential/${cohort}.diffbind_control_subtracted.log" 2>&1; then
                failures=$((failures+1)); mkdir -p "$root/sensitivity_control_subtracted"
                printf '{"status":"FAILED","module":"DiffBind_control_subtracted"}\n' > "$root/sensitivity_control_subtracted/FAILED.json"
            fi
        fi
    fi
    if is_true "$RUN_TARGET_CONTROL_INTERACTION"; then
        if ! run_logged Rscript "${SCRIPT_DIR}/control_interaction_analysis.R" "$SAMPLE_MANIFEST" "$cohort" "$consensus" "$OUTPUT_DIR" \
            "$root/sensitivity_target_control_interaction/deseq2" "$DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION" "$DIFFERENTIAL_ALPHA" \
            >"${OUTPUT_DIR}/logs/differential/${cohort}.control_interaction.log" 2>&1; then
            failures=$((failures+1)); mkdir -p "$root/sensitivity_target_control_interaction"
            printf '{"status":"FAILED","module":"target_control_interaction"}\n' > "$root/sensitivity_target_control_interaction/FAILED.json"
        fi
    fi
    printf 'analysis_family\trole\nprimary_target_only\tPRIMARY\nsensitivity_control_subtracted\tSENSITIVITY\nsensitivity_target_control_interaction\tSENSITIVITY\n' \
        > "$root/concordance/analysis_roles.tsv"
    printf '%s\t%s\t%s\n' "$cohort" "$failures" "$skips" > "$status_dir/${cohort}.tsv"
}

parallel_pool_init "$DIFFERENTIAL_PARALLEL_JOBS"
while IFS=$'\t' read -r cohort cohort_key genome assay_profile factor antibody_id layout target_class duplicate_policy caller peak_class n_samples sample_keys conditions; do
    [[ "$cohort" == "cohort_id" ]] && continue
    parallel_pool_submit "differential:$cohort" differential_worker "$cohort" "$caller" "$peak_class"
done < "$COHORT_MANIFEST"
parallel_pool_wait_all

failures="$(awk -F '\t' '{n+=$2} END{print n+0}' "$status_dir"/*.tsv 2>/dev/null || echo 0)"
skips="$(awk -F '\t' '{n+=$3} END{print n+0}' "$status_dir"/*.tsv 2>/dev/null || echo 0)"
if (( failures > 0 )); then differential_status=FAILED
elif (( skips > 0 )); then differential_status=COMPLETED_WITH_WARNINGS
else differential_status=SUCCESS
fi
printf 'status\tfailed_modules\tskipped_cohorts\tparallel_jobs\n%s\t%s\t%s\t%s\n' \
    "$differential_status" "$failures" "$skips" "$DIFFERENTIAL_PARALLEL_JOBS" > "${OUTPUT_DIR}/08_differential/stage_status.tsv"
(( failures == 0 ))

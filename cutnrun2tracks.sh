#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$(<"${SCRIPT_DIR}/VERSION")"
CONFIG=""
SAMPLESHEET_OVERRIDE=""
OUTPUT_OVERRIDE=""
PLAN_ONLY=false
PREFLIGHT_ONLY=false
FROM_STAGE=""
STOP_AFTER=""
SHOW_VERSION=false

usage() {
    cat <<'USAGE'
Usage: cutnrun2tracks.sh --config config.conf [options]

Options:
  --samplesheet FILE    Override SAMPLESHEET from config
  --output-dir DIR      Override OUTPUT_DIR from config
  --plan                Validate metadata and write a stage plan; do not require bioinformatics tools
  --preflight-only      Validate metadata, tools, and references, then stop
  --from-stage NAME     Reuse validated earlier outputs; re-run NAME and all later stages
  --stop-after NAME     Stop after the named stage
  --version             Print workflow version and stop
  -h, --help            Show this help

Stages: preflight preprocess alignment filtering cpm peakcalling consensus
        spikein normalized_tracks metagene qc differential annotation report cleanup finalize
USAGE
}

while (( $# )); do
    case "$1" in
        --config|-c) CONFIG="${2:?missing config path}"; shift 2 ;;
        --samplesheet|-s) SAMPLESHEET_OVERRIDE="${2:?missing samplesheet path}"; shift 2 ;;
        --output-dir|-o) OUTPUT_OVERRIDE="${2:?missing output path}"; shift 2 ;;
        --plan) PLAN_ONLY=true; shift ;;
        --preflight-only) PREFLIGHT_ONLY=true; shift ;;
        --from-stage) FROM_STAGE="${2:?missing stage}"; shift 2 ;;
        --stop-after) STOP_AFTER="${2:?missing stage}"; shift 2 ;;
        --version) SHOW_VERSION=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
is_true() { [[ "${1,,}" == "true" ]]; }
if is_true "$SHOW_VERSION"; then printf 'cutnrun2tracks %s\n' "$VERSION"; exit 0; fi
[[ -n "$CONFIG" && -f "$CONFIG" ]] || { echo "ERROR: --config must name a readable file" >&2; exit 2; }
(( BASH_VERSINFO[0] > 5 || (BASH_VERSINFO[0] == 5 && BASH_VERSINFO[1] >= 1) )) || {
    echo "ERROR: Bash >=5.1 is required" >&2; exit 2;
}

temporary="$(mktemp -d)"
trap 'rm -rf "$temporary"' EXIT
python3 "${SCRIPT_DIR}/scripts/sanitize_text_inputs.py" "$CONFIG"
config_args=("$CONFIG" --template "${SCRIPT_DIR}/config/config.conf.template")
[[ -n "$SAMPLESHEET_OVERRIDE" ]] && config_args+=(--samplesheet "$SAMPLESHEET_OVERRIDE")
[[ -n "$OUTPUT_OVERRIDE" ]] && config_args+=(--output-dir "$OUTPUT_OVERRIDE")
python3 "${SCRIPT_DIR}/scripts/validate_config.py" "${config_args[@]}" --write-shell "$temporary/resolved.conf"
# shellcheck disable=SC1090
source "$temporary/resolved.conf"
python3 "${SCRIPT_DIR}/scripts/sanitize_text_inputs.py" "$SAMPLESHEET"
mkdir -p "$OUTPUT_DIR/00_metadata" "$OUTPUT_DIR/.checkpoints" "$OUTPUT_DIR/logs"
python3 "${SCRIPT_DIR}/scripts/validate_config.py" "${config_args[@]}" \
    --write-shell "$OUTPUT_DIR/00_metadata/resolved_config.conf" --write-tsv "$OUTPUT_DIR/00_metadata/resolved_config.tsv"
export C2T_CONFIG="$OUTPUT_DIR/00_metadata/resolved_config.conf"
# shellcheck disable=SC1090
source "$C2T_CONFIG"
if is_true "$WRITE_CONSOLE_LOG" && [[ "${C2T_CONSOLE_CAPTURED:-false}" != "true" ]]; then
    export C2T_CONSOLE_CAPTURED=true
    exec > >(tee -a "$OUTPUT_DIR/logs/cutnrun2tracks.console.log") 2>&1
fi
cp "$SAMPLESHEET" "$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv"

sample_args=("$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv"
    --spikein-mode "$SPIKEIN_MODE" --spikein-reference-id "$SPIKEIN_REFERENCE_ID"
    --peak-callers "$PEAK_CALLERS" --primary-peak-caller "$PRIMARY_PEAK_CALLER"
    --output-dir "$OUTPUT_DIR/00_metadata")
while IFS= read -r genome; do
    [[ -n "$genome" ]] || continue
    reference_key="BLACKLIST_$(printf '%s' "$genome" | tr '[:lower:].-' '[:upper:]__')"
    reference_path="${!reference_key:-}"
    [[ -n "$reference_path" ]] || { echo "ERROR: missing reference setting $reference_key" >&2; exit 2; }
    sample_args+=(--blacklist-map "$genome=$reference_path")
done < <(python3 -c 'import csv,sys; print("\n".join(sorted({r["genome"].strip() for r in csv.DictReader(open(sys.argv[1], encoding="utf-8", newline="")) if r.get("genome", "").strip()})))' "$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv")
is_true "$ALLOW_SHARED_CONTROLS" && sample_args+=(--allow-shared-controls)
is_true "$ALLOW_MIXED_LAYOUTS" && sample_args+=(--allow-mixed-layouts)
is_true "$ALLOW_MIXED_GENOMES" && sample_args+=(--allow-mixed-genomes)
is_true "$ALLOW_CONTROL_FREE_PEAKCALL" && sample_args+=(--allow-control-free)
if ! is_true "$PLAN_ONLY"; then sample_args+=(--check-files); fi
python3 "${SCRIPT_DIR}/scripts/validate_samplesheet.py" "${sample_args[@]}"
RUN_ASSAY_PROFILES="$(awk -F '\t' 'NR>1 {seen[$6]=1} END {for (value in seen) print value}' "$OUTPUT_DIR/00_metadata/sample_manifest.tsv" | LC_ALL=C sort | paste -sd, -)"
[[ -n "$RUN_ASSAY_PROFILES" ]] || { echo "ERROR: no assay profiles resolved from samplesheet" >&2; exit 2; }

stages=(preflight preprocess alignment filtering cpm peakcalling consensus spikein normalized_tracks metagene qc differential annotation report cleanup finalize)
for requested_stage in "$FROM_STAGE" "$STOP_AFTER"; do
    [[ -z "$requested_stage" ]] && continue
    valid=false
    for stage_name in "${stages[@]}"; do [[ "$requested_stage" == "$stage_name" ]] && valid=true; done
    is_true "$valid" || { echo "ERROR: unknown stage: $requested_stage" >&2; exit 2; }
done
if is_true "$PLAN_ONLY"; then
    printf 'order\tstage\n' > "$OUTPUT_DIR/00_metadata/planned_stages.tsv"
    for index in "${!stages[@]}"; do printf '%s\t%s\n' "$((index+1))" "${stages[$index]}" >> "$OUTPUT_DIR/00_metadata/planned_stages.tsv"; done
    printf 'cutnrun2tracks %s plan validated: %s\n' "$VERSION" "$OUTPUT_DIR"
    exit 0
fi

python3 "${SCRIPT_DIR}/scripts/reference_manifest.py" "$C2T_CONFIG" \
    "$OUTPUT_DIR/00_metadata/sample_manifest.tsv" "$OUTPUT_DIR/00_metadata/reference_manifest.tsv"
RUN_SIGNATURE="$(python3 "${SCRIPT_DIR}/scripts/checkpoint.py" signature --jobs "$CHECKPOINT_PARALLEL_JOBS" \
    "$SCRIPT_DIR/VERSION" "$SCRIPT_DIR/scripts" "$SCRIPT_DIR/common" "$C2T_CONFIG" \
    "$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv" "$OUTPUT_DIR/00_metadata/reference_manifest.tsv")"
printf 'run_id\tworkflow_version\trun_signature\tassay_profiles\tspikein_mode\n%s\t%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "$VERSION" "$RUN_SIGNATURE" "$RUN_ASSAY_PROFILES" "$SPIKEIN_MODE" > "$OUTPUT_DIR/00_metadata/run_manifest.tsv"
if is_true "$WRITE_COMMAND_LOG"; then
    touch "$OUTPUT_DIR/00_metadata/commands.log"
    if [[ ! -s "$OUTPUT_DIR/00_metadata/command_events.tsv" ]]; then
        printf 'command_id\tstart_utc\tend_utc\telapsed_seconds\texit_code\tcommand\n' \
            > "$OUTPUT_DIR/00_metadata/command_events.tsv"
    fi
else
    rm -f -- "$OUTPUT_DIR/00_metadata/commands.log" "$OUTPUT_DIR/00_metadata/command_events.tsv"
fi

timing_table="$OUTPUT_DIR/00_metadata/stage_timing.tsv"
events_table="$OUTPUT_DIR/00_metadata/workflow_events.tsv"
if [[ ! -s "$timing_table" ]]; then
    printf 'run_id\tstage\tstatus\tstart_utc\tend_utc\telapsed_seconds\n' > "$timing_table"
fi
if [[ ! -s "$events_table" ]]; then
    printf 'timestamp_utc\trun_id\tstage\tscope\tevent\texit_code\telapsed_seconds\tmessage\n' > "$events_table"
fi
record_event() {
    is_true "$WRITE_STRUCTURED_LOG" || return 0
    local stage="$1" scope="$2" event="$3" exit_code="$4" elapsed="$5" message="$6"
    message="${message//$'\t'/ }"; message="${message//$'\n'/ }"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        "$RUN_ID" "$stage" "$scope" "$event" "$exit_code" "$elapsed" "$message" >> "$events_table"
}
record_stage_timing() {
    local stage="$1" status="$2" start_utc="$3" start_epoch="$4"
    local end_utc end_epoch elapsed
    end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; end_epoch="$(date -u +%s)"; elapsed=$((end_epoch-start_epoch))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$RUN_ID" "$stage" "$status" "$start_utc" "$end_utc" "$elapsed" >> "$timing_table"
    printf '%s' "$elapsed"
}
CURRENT_STAGE=initialization
PIPELINE_FINISHED=false
on_exit() {
    local status="$1"
    rm -rf -- "$temporary"
    if (( status != 0 )) && ! is_true "$PIPELINE_FINISHED"; then
        record_event "$CURRENT_STAGE" workflow FAILED "$status" 0 "workflow terminated unexpectedly"
        echo "=== [$CURRENT_STAGE] WORKFLOW FAILED === exit=$status" >&2
    fi
}
on_signal() {
    local signal="$1"
    record_event "$CURRENT_STAGE" workflow "SIGNAL_$signal" 130 0 "workflow interrupted"
    exit 130
}
trap 'on_exit $?' EXIT
trap 'on_signal INT' INT
trap 'on_signal TERM' TERM
trap 'on_signal HUP' HUP

force=false
before_from_stage=false
[[ -n "$FROM_STAGE" ]] && before_from_stage=true
run_stage() {
    local stage="$1"
    local command="$2"
    local checkpoint="$OUTPUT_DIR/.checkpoints/${stage}.json"
    shift 2
    local outputs=("$@") start_utc start_epoch elapsed status
    start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_epoch="$(date -u +%s)"
    CURRENT_STAGE="$stage"
    if [[ "$stage" == "$FROM_STAGE" ]]; then
        before_from_stage=false
        force=true
    fi
    if is_true "$before_from_stage"; then
        if python3 "$SCRIPT_DIR/scripts/checkpoint.py" adopt --checkpoint "$checkpoint" \
            --stage "$stage" --signature "$RUN_SIGNATURE" --jobs "$CHECKPOINT_PARALLEL_JOBS"; then
            echo "=== [$stage] REUSED: validated prior-stage outputs (--from-stage $FROM_STAGE) ==="
            status=REUSED
        else
            elapsed="$(record_stage_timing "$stage" FAILED_REUSE "$start_utc" "$start_epoch")"
            record_event "$stage" stage FAILED_REUSE 1 "$elapsed" "invalid or missing prior-stage checkpoint"
            echo "ERROR: cannot reuse invalid or missing prior-stage checkpoint: $stage" >&2
            exit 1
        fi
    elif [[ "$stage" != "cleanup" ]] && ! is_true "$force" && python3 "$SCRIPT_DIR/scripts/checkpoint.py" check --checkpoint "$checkpoint" --stage "$stage" --signature "$RUN_SIGNATURE" --jobs "$CHECKPOINT_PARALLEL_JOBS"; then
        echo "=== [$stage] SKIPPED: valid signature-and-output checkpoint ==="
        status=SKIPPED
    else
        echo "=== [$stage] START === $start_utc"
        record_event "$stage" stage START 0 0 "$command"
        if ! bash "$command"; then
            elapsed="$(record_stage_timing "$stage" FAILED "$start_utc" "$start_epoch")"
            record_event "$stage" stage FAILED 1 "$elapsed" "stage command failed"
            echo "=== [$stage] FAILED === elapsed=${elapsed}s" >&2
            return 1
        fi
        if ! python3 "$SCRIPT_DIR/scripts/checkpoint.py" write --checkpoint "$checkpoint" --stage "$stage" \
                --signature "$RUN_SIGNATURE" --outputs "${outputs[@]}" --jobs "$CHECKPOINT_PARALLEL_JOBS"; then
            elapsed="$(record_stage_timing "$stage" FAILED_CHECKPOINT "$start_utc" "$start_epoch")"
            record_event "$stage" stage FAILED_CHECKPOINT 1 "$elapsed" "checkpoint creation failed"
            return 1
        fi
        status=COMPLETE
    fi
    elapsed="$(record_stage_timing "$stage" "$status" "$start_utc" "$start_epoch")"
    record_event "$stage" stage "$status" 0 "$elapsed" "."
    [[ "$status" == "COMPLETE" ]] && echo "=== [$stage] COMPLETE === elapsed=${elapsed}s"
    if [[ "$stage" == "$STOP_AFTER" ]]; then echo "Stopped after $stage"; exit 0; fi
}

run_stage preflight "$SCRIPT_DIR/scripts/preflight.sh" "$OUTPUT_DIR/00_metadata/preflight_status.tsv"
is_true "$PREFLIGHT_ONLY" && { echo "Preflight complete"; exit 0; }
run_stage preprocess "$SCRIPT_DIR/scripts/preprocess_batch.sh" "$OUTPUT_DIR/02_trimmed_fastq"
run_stage alignment "$SCRIPT_DIR/scripts/align_batch.sh" "$OUTPUT_DIR/03_alignment/sorted"
run_stage filtering "$SCRIPT_DIR/scripts/mark_filter_batch.sh" "$OUTPUT_DIR/03_alignment/analysis"
run_stage cpm "$SCRIPT_DIR/scripts/coverage_batch.sh" "$OUTPUT_DIR/04_tracks/cpm"
run_stage peakcalling "$SCRIPT_DIR/scripts/peakcall_batch.sh" "$OUTPUT_DIR/05_peaks/per_sample"
run_stage consensus "$SCRIPT_DIR/scripts/consensus_batch.sh" "$OUTPUT_DIR/05_peaks/consensus"
run_stage spikein "$SCRIPT_DIR/scripts/spikein_batch.sh" "$OUTPUT_DIR/04_tracks/spikein"
run_stage normalized_tracks "$SCRIPT_DIR/scripts/normalized_tracks_batch.sh" "$OUTPUT_DIR/04_tracks"
run_stage metagene "$SCRIPT_DIR/scripts/metagene_batch.sh" "$OUTPUT_DIR/06_qc/metagene/status.tsv"
run_stage qc "$SCRIPT_DIR/scripts/qc_batch.sh" "$OUTPUT_DIR/06_qc"
run_stage differential "$SCRIPT_DIR/scripts/differential_batch.sh" "$OUTPUT_DIR/08_differential"
run_stage annotation "$SCRIPT_DIR/scripts/annotate_browser.sh" "$OUTPUT_DIR/07_annotation" "$OUTPUT_DIR/09_browser"
run_stage report "$SCRIPT_DIR/scripts/report_batch.sh" "$OUTPUT_DIR/10_reports"
run_stage cleanup "$SCRIPT_DIR/scripts/cleanup.sh" "$OUTPUT_DIR/00_metadata/cleanup_status.tsv" "$OUTPUT_DIR/00_metadata/cleanup_manifest.tsv"

CURRENT_STAGE=finalize
finalize_start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; finalize_start_epoch="$(date -u +%s)"
record_event finalize stage START 0 0 "final checksums"
if is_true "$WRITE_FILE_CHECKSUMS"; then
    python3 "$SCRIPT_DIR/scripts/finalize.py" "$OUTPUT_DIR" --jobs "$CHECKSUM_PARALLEL_JOBS"
else
    printf 'checksums disabled\n' > "$OUTPUT_DIR/00_metadata/final_checksums.sha256"
fi
python3 "$SCRIPT_DIR/scripts/checkpoint.py" write --checkpoint "$OUTPUT_DIR/.checkpoints/finalize.json" --stage finalize \
    --signature "$RUN_SIGNATURE" --outputs "$OUTPUT_DIR/00_metadata/final_checksums.sha256" --jobs "$CHECKPOINT_PARALLEL_JOBS"
finalize_elapsed="$(record_stage_timing finalize COMPLETE "$finalize_start_utc" "$finalize_start_epoch")"
record_event finalize stage COMPLETE 0 "$finalize_elapsed" "."
echo "=== [finalize] COMPLETE === elapsed=${finalize_elapsed}s"
if [[ "$STOP_AFTER" == "finalize" ]]; then echo "Stopped after finalize"; exit 0; fi
if grep -q '^COMPLETED_WITH_WARNINGS' "$OUTPUT_DIR/05_peaks/per_sample/stage_status.tsv" 2>/dev/null || \
    awk -F '\t' 'NR>1 && $2!="SUCCESS" {found=1} END{exit !found}' \
        "$OUTPUT_DIR/05_peaks/consensus/consensus_status.tsv" 2>/dev/null || \
    awk 'NR>1 {found=1} END{exit !found}' "$OUTPUT_DIR/10_reports/warning_summary.tsv" 2>/dev/null; then
    echo "cutnrun2tracks $VERSION completed with warnings: $OUTPUT_DIR"
else
    echo "cutnrun2tracks $VERSION completed successfully: $OUTPUT_DIR"
fi
PIPELINE_FINISHED=true

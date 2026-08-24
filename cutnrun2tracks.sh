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
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
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
cp "$SAMPLESHEET" "$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv"

sample_args=("$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv" --assay-profile "$ASSAY_PROFILE"
    --spikein-mode "$SPIKEIN_MODE" --spikein-reference-id "$SPIKEIN_REFERENCE_ID"
    --peak-callers "$PEAK_CALLERS" --primary-peak-caller "$PRIMARY_PEAK_CALLER"
    --output-dir "$OUTPUT_DIR/00_metadata")
is_true() { [[ "${1,,}" == "true" ]]; }
is_true "$ALLOW_SHARED_CONTROLS" && sample_args+=(--allow-shared-controls)
is_true "$ALLOW_MIXED_LAYOUTS" && sample_args+=(--allow-mixed-layouts)
is_true "$ALLOW_MIXED_GENOMES" && sample_args+=(--allow-mixed-genomes)
is_true "$ALLOW_CONTROL_FREE_PEAKCALL" && sample_args+=(--allow-control-free)
if ! is_true "$PLAN_ONLY"; then sample_args+=(--check-files); fi
python3 "${SCRIPT_DIR}/scripts/validate_samplesheet.py" "${sample_args[@]}"

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
RUN_SIGNATURE="$(python3 "${SCRIPT_DIR}/scripts/checkpoint.py" signature \
    "$SCRIPT_DIR/VERSION" "$SCRIPT_DIR/scripts" "$SCRIPT_DIR/common" "$C2T_CONFIG" \
    "$OUTPUT_DIR/00_metadata/sanitized_samplesheet.csv" "$OUTPUT_DIR/00_metadata/reference_manifest.tsv")"
printf 'run_id\tworkflow_version\trun_signature\tassay_profile\tspikein_mode\n%s\t%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "$VERSION" "$RUN_SIGNATURE" "$ASSAY_PROFILE" "$SPIKEIN_MODE" > "$OUTPUT_DIR/00_metadata/run_manifest.tsv"
: > "$OUTPUT_DIR/00_metadata/commands.log"

force=false
before_from_stage=false
[[ -n "$FROM_STAGE" ]] && before_from_stage=true
run_stage() {
    local stage="$1"
    local command="$2"
    local checkpoint="$OUTPUT_DIR/.checkpoints/${stage}.json"
    shift 2
    local outputs=("$@")
    if [[ "$stage" == "$FROM_STAGE" ]]; then
        before_from_stage=false
        force=true
    fi
    if is_true "$before_from_stage"; then
        if python3 "$SCRIPT_DIR/scripts/checkpoint.py" adopt --checkpoint "$checkpoint" \
            --stage "$stage" --signature "$RUN_SIGNATURE"; then
            echo "=== [$stage] REUSED: validated prior-stage outputs (--from-stage $FROM_STAGE) ==="
        else
            echo "ERROR: cannot reuse invalid or missing prior-stage checkpoint: $stage" >&2
            exit 1
        fi
    elif ! is_true "$force" && python3 "$SCRIPT_DIR/scripts/checkpoint.py" check --checkpoint "$checkpoint" --stage "$stage" --signature "$RUN_SIGNATURE"; then
        echo "=== [$stage] SKIPPED: valid signature-and-output checkpoint ==="
    else
        echo "=== [$stage] START ==="
        bash "$command"
        python3 "$SCRIPT_DIR/scripts/checkpoint.py" write --checkpoint "$checkpoint" --stage "$stage" \
            --signature "$RUN_SIGNATURE" --outputs "${outputs[@]}"
        echo "=== [$stage] COMPLETE ==="
    fi
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

if is_true "$WRITE_FILE_CHECKSUMS"; then
    python3 "$SCRIPT_DIR/scripts/finalize.py" "$OUTPUT_DIR"
else
    printf 'checksums disabled\n' > "$OUTPUT_DIR/00_metadata/final_checksums.sha256"
fi
python3 "$SCRIPT_DIR/scripts/checkpoint.py" write --checkpoint "$OUTPUT_DIR/.checkpoints/finalize.json" --stage finalize \
    --signature "$RUN_SIGNATURE" --outputs "$OUTPUT_DIR/00_metadata/final_checksums.sha256"
if [[ "$STOP_AFTER" == "finalize" ]]; then echo "Stopped after finalize"; exit 0; fi
echo "cutnrun2tracks $VERSION completed successfully: $OUTPUT_DIR"

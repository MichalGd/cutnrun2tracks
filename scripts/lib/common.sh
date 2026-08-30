#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
note() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

require_config() {
    [[ -n "${C2T_CONFIG:-}" && -f "$C2T_CONFIG" ]] || die "C2T_CONFIG is not a readable file"
    # shellcheck disable=SC1090
    source "$C2T_CONFIG"
    : "${OUTPUT_DIR:?OUTPUT_DIR missing}"
    SAMPLE_MANIFEST="${OUTPUT_DIR}/00_metadata/sample_manifest.tsv"
    COHORT_MANIFEST="${OUTPUT_DIR}/00_metadata/cohort_manifest.tsv"
    [[ -s "$SAMPLE_MANIFEST" ]] || die "sample manifest missing: $SAMPLE_MANIFEST"
}

is_true() { [[ "${1,,}" == "true" ]]; }

reference_value() {
    local prefix="${1:?prefix}" genome="${2:?genome}" key value
    key="${prefix}_$(printf '%s' "$genome" | tr '[:lower:].-' '[:upper:]__')"
    value="${!key:-}"
    [[ -n "$value" ]] || die "missing reference setting $key"
    printf '%s' "$value"
}

optional_reference_value() {
    local prefix="${1:?prefix}" genome="${2:?genome}" key
    key="${prefix}_$(printf '%s' "$genome" | tr '[:lower:].-' '[:upper:]__')"
    printf '%s' "${!key:-}"
}

signal_exclude_mask() {
    local layout="${1:?layout}" duplicate_policy="${2:-remove}"
    case "$layout:$duplicate_policy" in
        PE:retain) printf '2828' ;;
        PE:remove) printf '3852' ;;
        SE:retain) printf '2820' ;;
        SE:remove) printf '3844' ;;
        *) die "unsupported signal-count layout/policy: $layout/$duplicate_policy" ;;
    esac
}

signal_count() {
    local bam="${1:?bam}" layout="${2:?layout}" duplicate_policy="${3:-remove}" exclude
    exclude="$(signal_exclude_mask "$layout" "$duplicate_policy")"
    if [[ "$layout" == "PE" ]]; then
        samtools view -c -f 66 -F "$exclude" "$bam"
    else
        samtools view -c -F "$exclude" "$bam"
    fi
}

analysis_bam_path() {
    printf '%s/03_alignment/analysis/%s.host.analysis.bam' "$OUTPUT_DIR" "${1:?sample_key}"
}

policy_bam_path() {
    local key="${1:?sample_key}" policy="${2:?policy}"
    case "$policy" in
        analysis) analysis_bam_path "$key" ;;
        permissive) printf '%s/03_alignment/filtered/q0_dup-retained/%s.host.q0.dup-retained.bam' "$OUTPUT_DIR" "$key" ;;
        intermediate) printf '%s/03_alignment/filtered/q0_dup-removed/%s.host.q0.dup-removed.bam' "$OUTPUT_DIR" "$key" ;;
        stringent) printf '%s/03_alignment/filtered/q30_dup-removed/%s.host.q30.dup-removed.bam' "$OUTPUT_DIR" "$key" ;;
        *) die "unknown BAM policy: $policy" ;;
    esac
}

record_command() {
    is_true "${WRITE_COMMAND_LOG:-true}" || return 0
    local quoted=() argument
    for argument in "$@"; do quoted+=("$(printf '%q' "$argument")"); done
    printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${quoted[*]}" >> "${OUTPUT_DIR}/00_metadata/commands.log"
}

run_logged() {
    local quoted=() argument command_text command_id start_utc start_epoch end_utc end_epoch status elapsed
    for argument in "$@"; do quoted+=("$(printf '%q' "$argument")"); done
    command_text="${quoted[*]}"
    command_id="$(date -u +%s).${BASHPID}.${RANDOM}"
    start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_epoch="$(date -u +%s)"
    record_command "$@"
    if "$@"; then status=0; else status=$?; fi
    end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; end_epoch="$(date -u +%s)"; elapsed=$((end_epoch-start_epoch))
    if is_true "${WRITE_COMMAND_LOG:-true}"; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$command_id" "$start_utc" "$end_utc" "$elapsed" "$status" "$command_text" \
            >> "${OUTPUT_DIR}/00_metadata/command_events.tsv"
    fi
    return "$status"
}

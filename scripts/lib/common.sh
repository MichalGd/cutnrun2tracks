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

signal_count() {
    local bam="${1:?bam}" layout="${2:?layout}"
    if [[ "$layout" == "PE" ]]; then
        samtools view -c -f 66 -F 3840 "$bam"
    else
        samtools view -c -F 3844 "$bam"
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
    local quoted=() argument
    for argument in "$@"; do quoted+=("$(printf '%q' "$argument")"); done
    printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${quoted[*]}" >> "${OUTPUT_DIR}/00_metadata/commands.log"
}

run_logged() {
    record_command "$@"
    "$@"
}

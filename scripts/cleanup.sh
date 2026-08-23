#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

manifest="${OUTPUT_DIR}/00_metadata/cleanup_manifest.tsv"
printf 'path\taction\n' > "$manifest"
is_true "$ENABLE_AUTOMATIC_CLEANUP" || { printf 'status\nSKIPPED_DISABLED\n' > "${OUTPUT_DIR}/00_metadata/cleanup_status.tsv"; exit 0; }
[[ -s "${OUTPUT_DIR}/10_reports/pipeline_report.html" ]] || die "cleanup requires a validated final report"
resolved_output="$(realpath -m "$OUTPUT_DIR")"
[[ "$resolved_output" != "/" && "$resolved_output" != "." && ${#resolved_output} -gt 8 ]] || die "unsafe output root for cleanup: $resolved_output"
remove_child() {
    local target="$1" resolved
    [[ -e "$target" || -L "$target" ]] || return 0
    resolved="$(realpath -m "$target")"
    [[ "$resolved" == "$resolved_output"/* ]] || die "cleanup target escaped output root: $resolved"
    printf '%s\tdeleted_regenerable_intermediate\n' "$resolved" >> "$manifest"
    rm -rf -- "$resolved"
}

is_true "$KEEP_TRIMMED_FASTQ" || remove_child "${OUTPUT_DIR}/02_trimmed_fastq"
is_true "$KEEP_RAW_ALIGNMENT_BAMS" || remove_child "${OUTPUT_DIR}/03_alignment/sorted"
is_true "$KEEP_MARKED_BAMS" || remove_child "${OUTPUT_DIR}/03_alignment/marked"
if ! is_true "$KEEP_FILTERED_BAMS"; then
    # Preserve self-contained analysis BAMs before removing their symlink targets.
    find "${OUTPUT_DIR}/03_alignment/analysis" -type l -name '*.bam' -print | while read -r link; do
        target="$(readlink -f "$link")"; temporary="${link}.materialized"
        cp --reflink=auto "$target" "$temporary" 2>/dev/null || cp "$target" "$temporary"
        rm -f "$link" "${link}.bai"; mv "$temporary" "$link"
        cp "${target}.bai" "${link}.bai"
    done
    remove_child "${OUTPUT_DIR}/03_alignment/filtered"
fi
is_true "$KEEP_SPIKEIN_BAMS" || remove_child "${OUTPUT_DIR}/03_alignment/spikein/composite"
is_true "$KEEP_SPIKEIN_BAMS" || remove_child "${OUTPUT_DIR}/03_alignment/spikein/spike"
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/00_metadata/cleanup_status.tsv"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:?usage: generate_multiqc_report.sh OUTPUT_DIR [REPORT_DIR]}"
REPORT_DIR="${2:-${OUTPUT_DIR}/10_reports}"
REPORT_NAME=cutnrun2tracks_multiqc_report
REPORT_HTML="${REPORT_DIR}/${REPORT_NAME}.html"
REPORT_LOG="${REPORT_DIR}/multiqc.log"

[[ -d "$OUTPUT_DIR" ]] || { echo "ERROR: output directory missing: $OUTPUT_DIR" >&2; exit 2; }
[[ -s "${OUTPUT_DIR}/00_metadata/sample_manifest.tsv" ]] || {
    echo "ERROR: sample manifest missing under output directory: $OUTPUT_DIR" >&2
    exit 2
}
command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc is not available on PATH" >&2; exit 2; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is not available on PATH" >&2; exit 2; }

mkdir -p "$REPORT_DIR"
CUSTOM_DIR="$(mktemp -d "${TMPDIR:-/tmp}/cutnrun2tracks-multiqc.XXXXXX")"
cleanup_custom_dir() { rm -rf -- "$CUSTOM_DIR"; }
trap cleanup_custom_dir EXIT

python3 "${SCRIPT_DIR}/prepare_multiqc_content.py" "$OUTPUT_DIR" "$CUSTOM_DIR"
cp -- "${CUSTOM_DIR}/custom_content_manifest.tsv" \
    "${REPORT_DIR}/multiqc_custom_content_manifest.tsv"

# MultiQC 1.35 can reject or merge some native deepTools outputs. The workflow's
# authoritative deepTools files remain untouched; selected images and metagene
# tables are supplied through deterministic custom content instead.
if ! multiqc "$OUTPUT_DIR" "$CUSTOM_DIR" \
    --outdir "$REPORT_DIR" \
    --filename "$REPORT_NAME" \
    --exclude deeptools \
    --ignore "*/10_reports/*" \
    --cl-config 'ignore_images: false' \
    --data-format tsv --export --force \
    2>&1 | tee "$REPORT_LOG"; then
    echo "ERROR: MultiQC failed; inspect $REPORT_LOG" >&2
    exit 1
fi

if grep -Eiq 'Oops! The .* MultiQC module broke|ValidationError|Error converting colou?r' "$REPORT_LOG"; then
    echo "ERROR: MultiQC reported a module or validation failure; inspect $REPORT_LOG" >&2
    exit 1
fi
[[ -s "$REPORT_HTML" ]] || {
    echo "ERROR: MultiQC did not create a non-empty report: $REPORT_HTML" >&2
    exit 1
}

data_dir=""
for candidate in \
    "${REPORT_DIR}/${REPORT_NAME}_data" \
    "${REPORT_DIR}/multiqc_data"; do
    if [[ -d "$candidate" ]]; then data_dir="$candidate"; break; fi
done
[[ -n "$data_dir" ]] || {
    echo "ERROR: MultiQC report data directory was not created under $REPORT_DIR" >&2
    exit 1
}

custom_items="$(wc -l < "${REPORT_DIR}/multiqc_custom_content_manifest.tsv")"
custom_items=$(( custom_items > 0 ? custom_items - 1 : 0 ))
printf 'status\treport\tdata_dir\tcustom_content_items\nSUCCESS\t%s\t%s\t%s\n' \
    "$REPORT_HTML" "$data_dir" "$custom_items" > "${REPORT_DIR}/multiqc_status.tsv"
echo "Unified MultiQC report: $REPORT_HTML"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
    echo "Usage: regenerate_reports.sh --output-dir OUTPUT_DIR" >&2
}

OUTPUT_DIR=""
while (( $# )); do
    case "$1" in
        --output-dir) OUTPUT_DIR="${2:?missing output directory}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done
[[ -n "$OUTPUT_DIR" ]] || { usage; exit 2; }
[[ -s "${OUTPUT_DIR}/00_metadata/sample_manifest.tsv" ]] || {
    echo "ERROR: not a cutnrun2tracks output directory: $OUTPUT_DIR" >&2
    exit 2
}

python3 "${ROOT}/scripts/generate_report.py" "$OUTPUT_DIR"
bash "${ROOT}/scripts/generate_multiqc_report.sh" "$OUTPUT_DIR" "${OUTPUT_DIR}/10_reports"

(
    cd "${OUTPUT_DIR}/10_reports"
    find . -type f ! -name report_checksums.sha256 -print0 | sort -z |
        xargs -0 -r sha256sum > report_checksums.sha256
)

for required in \
    pipeline_report.html run_summary.tsv warning_summary.tsv \
    cutnrun2tracks_multiqc_report.html multiqc_status.tsv report_checksums.sha256; do
    [[ -s "${OUTPUT_DIR}/10_reports/${required}" ]] || {
        echo "ERROR: regenerated report artifact missing: $required" >&2
        exit 1
    }
done

echo "Reports regenerated without rerunning upstream analysis: ${OUTPUT_DIR}/10_reports"
echo "NOTE: 00_metadata/final_checksums.sha256 predates these regenerated reports; report_checksums.sha256 covers the new report bundle."

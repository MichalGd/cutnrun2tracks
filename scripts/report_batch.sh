#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config
run_logged python3 "${SCRIPT_DIR}/generate_report.py" "$OUTPUT_DIR"
if is_true "$RUN_MULTIQC"; then
    run_logged bash "${SCRIPT_DIR}/generate_multiqc_report.sh" "$OUTPUT_DIR" "${OUTPUT_DIR}/10_reports"
else
    printf 'status\treason\nSKIPPED\tRUN_MULTIQC=false\n' > "${OUTPUT_DIR}/10_reports/multiqc_status.tsv"
fi

for required in pipeline_report.html run_summary.tsv warning_summary.tsv; do
    [[ -s "${OUTPUT_DIR}/10_reports/${required}" ]] || \
        die "report artifact missing or empty: ${OUTPUT_DIR}/10_reports/${required}"
done
if is_true "$RUN_MULTIQC"; then
    [[ -s "${OUTPUT_DIR}/10_reports/cutnrun2tracks_multiqc_report.html" ]] || \
        die "unified MultiQC report missing"
fi

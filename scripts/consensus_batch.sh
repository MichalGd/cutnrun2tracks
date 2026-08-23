#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

mkdir -p "${OUTPUT_DIR}/05_peaks/consensus"
args=(
    --sample-manifest "$SAMPLE_MANIFEST"
    --cohort-manifest "$COHORT_MANIFEST"
    --output-root "${OUTPUT_DIR}/05_peaks/consensus"
    --minimum-support "$CONSENSUS_MIN_BIOLOGICAL_SAMPLES"
)
is_true "$ALLOW_SINGLE_SAMPLE_CONSENSUS" && args+=(--allow-single)
is_true "$REQUIRE_ALL_ENABLED_TRACKS" && args+=(--require-all)
run_logged python3 "${SCRIPT_DIR}/build_consensus.py" "${args[@]}"

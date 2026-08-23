#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

root="${OUTPUT_DIR}/06_qc/metagene"
mkdir -p "$root" "${OUTPUT_DIR}/00_metadata"
if ! is_true "$RUN_METAGENE"; then
    printf '{"status":"SKIPPED","reason":"RUN_METAGENE=false"}\n' > "${root}/SKIPPED.json"
    printf 'status\ttasks\nSKIPPED\t0\n' > "${root}/status.tsv"
    exit 0
fi

track_manifest="${OUTPUT_DIR}/00_metadata/metagene_tracks.tsv"
track_args=(
    --sample-manifest "$SAMPLE_MANIFEST"
    --resolved-config "$C2T_CONFIG"
    --output-dir "$OUTPUT_DIR"
    --track-family "$METAGENE_TRACK_FAMILY"
    --output "$track_manifest"
)
is_true "$METAGENE_ALLOW_CPM_FALLBACK" && track_args+=(--allow-cpm-fallback)
is_true "$METAGENE_INCLUDE_CONTROLS" && track_args+=(--include-controls)
run_logged python3 "${SCRIPT_DIR}/build_metagene_track_manifest.py" "${track_args[@]}"

runner=(
    bash "${REPOSITORY_ROOT}/common/metagene/run_metagene.sh"
    --track-manifest "$track_manifest"
    --gene-set-manifest "$METAGENE_GENE_SET_MANIFEST"
    --output-dir "$root"
    --gene-sets "$METAGENE_GENE_SETS"
    --modes "$METAGENE_MODES"
    --parallel-jobs "$METAGENE_PARALLEL_JOBS"
    --threads-per-job "$METAGENE_THREADS_COMPUTEMATRIX"
    --reference-upstream "$METAGENE_REFERENCE_UPSTREAM_BP"
    --reference-downstream "$METAGENE_REFERENCE_DOWNSTREAM_BP"
    --body-upstream "$METAGENE_BODY_UPSTREAM_BP"
    --body-downstream "$METAGENE_BODY_DOWNSTREAM_BP"
    --body-length "$METAGENE_BODY_LENGTH_BP"
    --bin-size "$METAGENE_BIN_SIZE_BP"
    --missing-data "$METAGENE_MISSING_DATA_POLICY"
    --color-map "$METAGENE_COLOR_MAP"
    --zmin "$METAGENE_ZMIN"
    --zmax "$METAGENE_ZMAX"
    --dpi "$METAGENE_DPI"
    --plot-formats "$METAGENE_PLOT_FORMATS"
)
is_true "$METAGENE_SKIP_ZERO_REGIONS" && runner+=(--skip-zero-regions)
run_logged "${runner[@]}"

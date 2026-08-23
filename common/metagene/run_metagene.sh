#!/usr/bin/env bash
set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "${MODULE_DIR}/../.." && pwd)"
# shellcheck source=../../scripts/lib/parallel_jobs.sh
source "${REPOSITORY_ROOT}/scripts/lib/parallel_jobs.sh"

TRACK_MANIFEST=""; GENE_SET_MANIFEST=""; OUTPUT_DIR=""
GENE_SETS="protein_coding,broadly_expressed"; MODES="tss,tes,gene_body"
PARALLEL_JOBS=2; THREADS_PER_JOB=4
REFERENCE_UPSTREAM=3000; REFERENCE_DOWNSTREAM=3000
BODY_UPSTREAM=3000; BODY_DOWNSTREAM=3000; BODY_LENGTH=5000; BIN_SIZE=10
MISSING_POLICY=zero; SKIP_ZERO=false; COLOR_MAP=viridis; ZMIN=auto; ZMAX=auto
DPI=200; PLOT_FORMATS=png,pdf; SKIP_BIGWIG_VALIDATION=false

usage() {
    cat <<'USAGE'
Usage: run_metagene.sh --track-manifest FILE --gene-set-manifest FILE --output-dir DIR [options]

Options:
  --gene-sets LIST              Comma-separated gene-set IDs
  --modes LIST                  Comma-separated tss,tes,gene_body
  --parallel-jobs N             Concurrent matrix jobs (default: 2)
  --threads-per-job N           computeMatrix threads per job (default: 4)
  --reference-upstream BP       TSS/TES upstream flank (default: 3000)
  --reference-downstream BP     TSS/TES downstream flank (default: 3000)
  --body-upstream BP            Scaled-body upstream flank (default: 3000)
  --body-downstream BP          Scaled-body downstream flank (default: 3000)
  --body-length BP              Scaled body length (default: 5000)
  --bin-size BP                 Matrix bin size (default: 10)
  --missing-data zero|na        Missing-value interpretation (default: zero)
  --skip-zero-regions           Exclude all-zero rows (not recommended)
  --color-map NAME              Matplotlib color map (default: viridis)
  --zmin VALUE                  Heatmap minimum or auto
  --zmax VALUE                  Heatmap maximum or auto
  --dpi N                       Raster plot DPI (default: 200)
  --plot-formats png,pdf        Required plot formats
  --skip-bigwig-validation      Testing/dry-run escape hatch
  -h, --help                    Show this help
USAGE
}

while (( $# )); do
    case "$1" in
        --track-manifest) TRACK_MANIFEST="${2:?missing value}"; shift 2 ;;
        --gene-set-manifest) GENE_SET_MANIFEST="${2:?missing value}"; shift 2 ;;
        --output-dir) OUTPUT_DIR="${2:?missing value}"; shift 2 ;;
        --gene-sets) GENE_SETS="${2:?missing value}"; shift 2 ;;
        --modes) MODES="${2:?missing value}"; shift 2 ;;
        --parallel-jobs) PARALLEL_JOBS="${2:?missing value}"; shift 2 ;;
        --threads-per-job) THREADS_PER_JOB="${2:?missing value}"; shift 2 ;;
        --reference-upstream) REFERENCE_UPSTREAM="${2:?missing value}"; shift 2 ;;
        --reference-downstream) REFERENCE_DOWNSTREAM="${2:?missing value}"; shift 2 ;;
        --body-upstream) BODY_UPSTREAM="${2:?missing value}"; shift 2 ;;
        --body-downstream) BODY_DOWNSTREAM="${2:?missing value}"; shift 2 ;;
        --body-length) BODY_LENGTH="${2:?missing value}"; shift 2 ;;
        --bin-size) BIN_SIZE="${2:?missing value}"; shift 2 ;;
        --missing-data) MISSING_POLICY="${2:?missing value}"; shift 2 ;;
        --skip-zero-regions) SKIP_ZERO=true; shift ;;
        --color-map) COLOR_MAP="${2:?missing value}"; shift 2 ;;
        --zmin) ZMIN="${2:?missing value}"; shift 2 ;;
        --zmax) ZMAX="${2:?missing value}"; shift 2 ;;
        --dpi) DPI="${2:?missing value}"; shift 2 ;;
        --plot-formats) PLOT_FORMATS="${2:?missing value}"; shift 2 ;;
        --skip-bigwig-validation) SKIP_BIGWIG_VALIDATION=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$TRACK_MANIFEST" && -f "$TRACK_MANIFEST" ]] || { echo "ERROR: readable --track-manifest is required" >&2; exit 2; }
[[ -n "$GENE_SET_MANIFEST" && -f "$GENE_SET_MANIFEST" ]] || { echo "ERROR: readable --gene-set-manifest is required" >&2; exit 2; }
[[ -n "$OUTPUT_DIR" ]] || { echo "ERROR: --output-dir is required" >&2; exit 2; }
for value in "$PARALLEL_JOBS" "$THREADS_PER_JOB" "$REFERENCE_UPSTREAM" "$REFERENCE_DOWNSTREAM" \
    "$BODY_UPSTREAM" "$BODY_DOWNSTREAM" "$BODY_LENGTH" "$BIN_SIZE" "$DPI"; do
    [[ "$value" =~ ^[0-9]+$ ]] || { echo "ERROR: numeric settings must be non-negative integers" >&2; exit 2; }
done
(( PARALLEL_JOBS > 0 && THREADS_PER_JOB > 0 && BODY_LENGTH > 0 && BIN_SIZE > 0 && DPI > 0 )) || {
    echo "ERROR: jobs, threads, body length, bin size, and DPI must be positive" >&2; exit 2;
}
[[ "$MISSING_POLICY" == "zero" || "$MISSING_POLICY" == "na" ]] || { echo "ERROR: --missing-data must be zero or na" >&2; exit 2; }
IFS=',' read -r -a requested_formats <<< "$PLOT_FORMATS"
format_set=","; for format in "${requested_formats[@]}"; do format_set+="${format,,},"; done
[[ "$format_set" == *",png,"* && "$format_set" == *",pdf,"* ]] || { echo "ERROR: --plot-formats must include png and pdf" >&2; exit 2; }
for format in "${requested_formats[@]}"; do
    [[ "${format,,}" == "png" || "${format,,}" == "pdf" ]] || { echo "ERROR: unsupported plot format: $format" >&2; exit 2; }
done
for command_name in python3 computeMatrix plotProfile plotHeatmap gzip; do
    command -v "$command_name" >/dev/null 2>&1 || { echo "ERROR: missing command: $command_name" >&2; exit 1; }
done

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/multiqc_data"
tasks="${OUTPUT_DIR}/tasks.tsv"
warnings="${OUTPUT_DIR}/input_warnings.tsv"
validator=(python3 "${MODULE_DIR}/validate_metagene_inputs.py"
    --track-manifest "$TRACK_MANIFEST" --gene-set-manifest "$GENE_SET_MANIFEST"
    --gene-sets "$GENE_SETS" --modes "$MODES" --output "$tasks" --warnings "$warnings")
[[ "$SKIP_BIGWIG_VALIDATION" == "true" ]] && validator+=(--skip-bigwig-validation)
"${validator[@]}"

export METAGENE_REFERENCE_UPSTREAM_BP="$REFERENCE_UPSTREAM"
export METAGENE_REFERENCE_DOWNSTREAM_BP="$REFERENCE_DOWNSTREAM"
export METAGENE_BODY_UPSTREAM_BP="$BODY_UPSTREAM"
export METAGENE_BODY_DOWNSTREAM_BP="$BODY_DOWNSTREAM"
export METAGENE_BODY_LENGTH_BP="$BODY_LENGTH"
export METAGENE_BIN_SIZE_BP="$BIN_SIZE"
export METAGENE_MISSING_DATA_POLICY="$MISSING_POLICY"
export METAGENE_SKIP_ZERO_REGIONS="$SKIP_ZERO"
export METAGENE_COLOR_MAP="$COLOR_MAP"
export METAGENE_ZMIN="$ZMIN"
export METAGENE_ZMAX="$ZMAX"
export METAGENE_DPI="$DPI"
export METAGENE_PLOT_FORMATS="$PLOT_FORMATS"
export METAGENE_THREADS_COMPUTEMATRIX="$THREADS_PER_JOB"

parallel_pool_init "$PARALLEL_JOBS"
while IFS=$'\t' read -r task_id sample_id assay genome bigwig normalization normalization_detail blacklist chrom_sizes \
    gene_set_id gene_set_label bed12 bed12_sha256 n_genes_reference mode; do
    [[ "$task_id" == "task_id" ]] && continue
    parallel_pool_submit "$task_id" bash "${MODULE_DIR}/metagene_worker.sh" \
        "$task_id" "$sample_id" "$assay" "$genome" "$bigwig" "$normalization" "$normalization_detail" \
        "$blacklist" "$gene_set_id" "$gene_set_label" "$bed12" "$bed12_sha256" "$n_genes_reference" "$mode" "$OUTPUT_DIR"
done < "$tasks"
parallel_pool_wait_all

python3 "${MODULE_DIR}/summarize_metagene.py" "$OUTPUT_DIR"
printf 'status\ttasks\nSUCCESS\t%s\n' "$(($(wc -l < "$tasks") - 1))" > "${OUTPUT_DIR}/status.tsv"
echo "Metagene analysis complete: $OUTPUT_DIR"

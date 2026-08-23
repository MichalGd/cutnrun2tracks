#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: metagene_worker.sh TASK_ID SAMPLE ASSAY GENOME BIGWIG NORMALIZATION NORMALIZATION_DETAIL BLACKLIST GENE_SET_ID GENE_SET_LABEL BED12 BED12_SHA256 N_GENES MODE OUTPUT_DIR" >&2
}

(( $# == 15 )) || { usage; exit 2; }
TASK_ID="$1"; SAMPLE_ID="$2"; ASSAY="$3"; GENOME="$4"; BIGWIG="$5"
NORMALIZATION="$6"; NORMALIZATION_DETAIL="$7"; BLACKLIST="$8"; GENE_SET_ID="$9"
GENE_SET_LABEL="${10}"; BED12="${11}"; BED12_SHA256="${12}"; N_GENES_REFERENCE="${13}"
MODE="${14}"; OUTPUT_DIR="${15}"

for required in METAGENE_REFERENCE_UPSTREAM_BP METAGENE_REFERENCE_DOWNSTREAM_BP \
    METAGENE_BODY_UPSTREAM_BP METAGENE_BODY_DOWNSTREAM_BP METAGENE_BODY_LENGTH_BP \
    METAGENE_BIN_SIZE_BP METAGENE_MISSING_DATA_POLICY METAGENE_SKIP_ZERO_REGIONS \
    METAGENE_COLOR_MAP METAGENE_ZMIN METAGENE_ZMAX METAGENE_DPI \
    METAGENE_PLOT_FORMATS METAGENE_THREADS_COMPUTEMATRIX; do
    [[ -n "${!required:-}" ]] || { echo "ERROR: missing worker setting $required" >&2; exit 2; }
done
[[ "$TASK_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "ERROR: unsafe task ID" >&2; exit 2; }
[[ -s "$BIGWIG" && -s "$BED12" && -s "$BLACKLIST" ]] || { echo "ERROR: task input is missing" >&2; exit 2; }
[[ "$MODE" == "tss" || "$MODE" == "tes" || "$MODE" == "gene_body" ]] || { echo "ERROR: invalid mode $MODE" >&2; exit 2; }

export MPLBACKEND=Agg
final_dir="${OUTPUT_DIR}/${GENOME}/${GENE_SET_ID}/${MODE}/${SAMPLE_ID}"
mkdir -p "$final_dir" "${OUTPUT_DIR}/logs"
temporary="$(mktemp -d "${final_dir}/.metagene.${TASK_ID}.XXXXXX")"
cleanup() { rm -rf -- "$temporary"; }
trap cleanup EXIT
prefix="${temporary}/${TASK_ID}"
command_log="${OUTPUT_DIR}/logs/${TASK_ID}.commands.log"
: > "$command_log"

run_cmd() {
    local quoted=() argument
    for argument in "$@"; do quoted+=("$(printf '%q' "$argument")"); done
    printf '%s\n' "${quoted[*]}" >> "$command_log"
    "$@"
}

optional_compute=()
[[ "${METAGENE_MISSING_DATA_POLICY,,}" == "zero" ]] && optional_compute+=(--missingDataAsZero)
[[ "${METAGENE_SKIP_ZERO_REGIONS,,}" == "true" ]] && optional_compute+=(--skipZeros)

common_compute=(
    --regionsFileName "$BED12"
    --scoreFileName "$BIGWIG"
    --binSize "$METAGENE_BIN_SIZE_BP"
    --averageTypeBins mean
    --blackListFileName "$BLACKLIST"
    --metagene
    --sortRegions keep
    --samplesLabel "$SAMPLE_ID"
    --outFileName "${prefix}.matrix.gz"
    --outFileNameMatrix "${prefix}.matrix.tsv"
    --outFileSortedRegions "${prefix}.matrix_regions.bed"
    --numberOfProcessors "$METAGENE_THREADS_COMPUTEMATRIX"
)

if [[ "$MODE" == "gene_body" ]]; then
    run_cmd computeMatrix scale-regions \
        "${common_compute[@]}" "${optional_compute[@]}" \
        --beforeRegionStartLength "$METAGENE_BODY_UPSTREAM_BP" \
        --regionBodyLength "$METAGENE_BODY_LENGTH_BP" \
        --afterRegionStartLength "$METAGENE_BODY_DOWNSTREAM_BP" \
        --startLabel TSS --endLabel TES
    mode_plot=(--startLabel TSS --endLabel TES)
    title_mode="scaled gene body"
else
    [[ "$MODE" == "tss" ]] && reference_point=TSS || reference_point=TES
    run_cmd computeMatrix reference-point \
        "${common_compute[@]}" "${optional_compute[@]}" \
        --referencePoint "$reference_point" \
        --beforeRegionStartLength "$METAGENE_REFERENCE_UPSTREAM_BP" \
        --afterRegionStartLength "$METAGENE_REFERENCE_DOWNSTREAM_BP"
    mode_plot=(--refPointLabel "$reference_point")
    title_mode="$reference_point"
fi

y_label="Signal (${NORMALIZATION}; normalized upstream)"
first_format=true
IFS=',' read -r -a formats <<< "$METAGENE_PLOT_FORMATS"
for format in "${formats[@]}"; do
    format="${format,,}"
    profile_extra=()
    heatmap_extra=()
    if [[ "$first_format" == "true" ]]; then
        profile_extra+=(--outFileNameData "${prefix}.profile.tsv")
        heatmap_extra+=(--outFileSortedRegions "${prefix}.heatmap_sorted_regions.bed")
        first_format=false
    fi
    run_cmd plotProfile \
        --matrixFile "${prefix}.matrix.gz" \
        --outFileName "${prefix}.profile.${format}" \
        --plotFileFormat "$format" \
        --plotType lines --averageType mean \
        "${mode_plot[@]}" \
        --regionsLabel "$GENE_SET_LABEL" --samplesLabel "$SAMPLE_ID" \
        --yAxisLabel "$y_label" \
        --plotTitle "$SAMPLE_ID | $GENE_SET_LABEL | $title_mode" \
        --dpi "$METAGENE_DPI" "${profile_extra[@]}"
    run_cmd plotHeatmap \
        --matrixFile "${prefix}.matrix.gz" \
        --outFileName "${prefix}.heatmap.${format}" \
        --plotFileFormat "$format" \
        --sortRegions descend --sortUsing mean \
        --colorMap "$METAGENE_COLOR_MAP" \
        --zMin "$METAGENE_ZMIN" --zMax "$METAGENE_ZMAX" \
        --whatToShow "heatmap and colorbar" \
        "${mode_plot[@]}" \
        --regionsLabel "$GENE_SET_LABEL" --samplesLabel "$SAMPLE_ID" \
        --plotTitle "$SAMPLE_ID | $GENE_SET_LABEL | $title_mode" \
        --dpi "$METAGENE_DPI" "${heatmap_extra[@]}"
done

gzip -f "${prefix}.matrix.tsv"
n_genes_matrix="$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "${prefix}.matrix_regions.bed")"
for output in "${prefix}".*; do
    mv -f -- "$output" "$final_dir/"
done

final_prefix="${final_dir}/${TASK_ID}"
metadata="${final_prefix}.metadata.tsv"
printf 'task_id\tsample_id\tassay\tgenome\tnormalization\tnormalization_detail\tgene_set\tmode\tn_genes_reference\tn_genes_matrix\tbigwig\tbed12\tbed12_sha256\tmatrix\tprofile_data\tprofile_png\tprofile_pdf\theatmap_png\theatmap_pdf\tsorted_regions\tstatus\n' > "$metadata"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS\n' \
    "$TASK_ID" "$SAMPLE_ID" "$ASSAY" "$GENOME" "$NORMALIZATION" "$NORMALIZATION_DETAIL" \
    "$GENE_SET_ID" "$MODE" "$N_GENES_REFERENCE" "$n_genes_matrix" "$BIGWIG" "$BED12" "$BED12_SHA256" \
    "${final_prefix}.matrix.gz" "${final_prefix}.profile.tsv" \
    "${final_prefix}.profile.png" "${final_prefix}.profile.pdf" \
    "${final_prefix}.heatmap.png" "${final_prefix}.heatmap.pdf" \
    "${final_prefix}.heatmap_sorted_regions.bed" >> "$metadata"
printf 'SUCCESS\n' > "${final_prefix}.SUCCESS"

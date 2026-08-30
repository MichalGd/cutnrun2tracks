#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/07_annotation" "${OUTPUT_DIR}/09_browser/ucsc" "${OUTPUT_DIR}/09_browser/igv"
trackdb="${OUTPUT_DIR}/09_browser/ucsc/trackDb.txt"
: > "$trackdb"

stream_text() {
    local path="${1:?path}"
    if [[ "$path" == *.gz ]]; then gzip -cd -- "$path"; else cat -- "$path"; fi
}

annotate_consensus() {
    local cohort="$1" genome="$2" caller="$3" peak_class="$4" consensus_dir consensus out gtf chrom_sizes genes genes_unsorted sorted_consensus ccre
    consensus_dir="${OUTPUT_DIR}/05_peaks/consensus/${cohort}/${caller}/${peak_class}"
    consensus="$(find "$consensus_dir" -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
    [[ -n "$consensus" ]] || return 0
    out="${OUTPUT_DIR}/07_annotation/${cohort}/consensus"; mkdir -p "$out"
    if is_true "$RUN_SIMPLE_PEAK_ANNOTATION"; then
        gtf="$(reference_value GTF "$genome")"; chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
        genes="${out}/genes.bed"; genes_unsorted="${out}/genes.unsorted.bed"
        sorted_consensus="${out}/${cohort}.consensus.sorted.bed"
        stream_text "$gtf" | awk 'BEGIN{FS=OFS="\t"} $3=="gene" {id="."; name="."; if(match($9,/gene_id "[^"]+"/))id=substr($9,RSTART+9,RLENGTH-10); if(match($9,/gene_name "[^"]+"/))name=substr($9,RSTART+11,RLENGTH-12); print $1,$4-1,$5,name,id,$7}' |
            awk 'NF==6 && $2>=0 && $3>$2' > "$genes_unsorted"
        bedtools sort -faidx "$chrom_sizes" -i "$genes_unsorted" > "$genes"
        bedtools sort -faidx "$chrom_sizes" -i "$consensus" > "$sorted_consensus"
        rm -f -- "$genes_unsorted"
        bedtools closest -a "$sorted_consensus" -b "$genes" -d -g "$chrom_sizes" > "${out}/${cohort}.nearest_gene.tsv"
    fi
    if is_true "$RUN_CCRE_ANNOTATION"; then
        ccre="$(optional_reference_value CCRE_BED "$genome")"
        if [[ -n "$ccre" && "$ccre" != "." && -s "$ccre" ]]; then
            bedtools intersect -wao -a "$consensus" -b "$ccre" > "${out}/${cohort}.ccre_reference_overlaps.tsv"
        else
            printf '{"status":"SKIPPED","reason":"CCRE reference unavailable"}\n' > "${out}/ccre_SKIPPED.json"
        fi
    fi
}

parallel_pool_init "$ANNOTATION_PARALLEL_JOBS"
while IFS=$'\t' read -r cohort cohort_key genome assay_profile factor antibody_id layout target_class duplicate_policy caller peak_class n_samples sample_keys conditions; do
    [[ "$cohort" == "cohort_id" ]] && continue
    parallel_pool_submit "annotation:$cohort" annotate_consensus "$cohort" "$genome" "$caller" "$peak_class"
done < "$COHORT_MANIFEST"
parallel_pool_wait_all

if is_true "$RUN_FEATURE_ANNOTATION_SUMMARY"; then
    feature_args=(--output-dir "$OUTPUT_DIR" --annotation-dir "$OUTPUT_DIR/07_annotation/feature_summary"
        --promoter-upstream "$PEAK_ANNOTATION_PROMOTER_UPSTREAM"
        --promoter-downstream "$PEAK_ANNOTATION_PROMOTER_DOWNSTREAM"
        --gene-end-window "$PEAK_ANNOTATION_GENE_END_WINDOW"
        --precedence "$PEAK_ANNOTATION_FEATURE_PRECEDENCE"
        --plot-formats "$PEAK_ANNOTATION_PLOT_FORMATS")
    is_true "$PEAK_ANNOTATION_INCLUDE_CONSENSUS" && feature_args+=(--include-consensus)
    while IFS= read -r genome; do
        [[ -n "$genome" ]] || continue
        feature_args+=(--gtf "$genome=$(reference_value GTF "$genome")")
        feature_args+=(--chrom-sizes "$genome=$(reference_value CHROM_SIZES "$genome")")
        ccre="$(optional_reference_value CCRE_BED "$genome")"
        [[ -n "$ccre" && "$ccre" != "." && -s "$ccre" ]] && feature_args+=(--ccre "$genome=$ccre")
    done < <(awk -F '\t' 'NR>1 {print $5}' "$SAMPLE_MANIFEST" | LC_ALL=C sort -u)
    run_logged python3 "${SCRIPT_DIR}/summarize_peak_annotations.py" "${feature_args[@]}"
fi

while IFS=$'\t' read -r sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition rest; do
    [[ "$sample_key" == "sample_key" ]] && continue
    bw="${OUTPUT_DIR}/04_tracks/cpm/${sample_key}.CPM.bw"
    [[ -s "$bw" ]] || continue
    url="$bw"
    [[ -n "$UCSC_BIGDATA_URL_BASE" ]] && url="${UCSC_BIGDATA_URL_BASE%/}/04_tracks/cpm/${sample_key}.CPM.bw"
    printf 'track %s_%s\ntype bigWig\nbigDataUrl %s\nshortLabel %s CPM\nlongLabel %s %s %s CPM\nvisibility full\nautoScale on\n\n' \
        "$UCSC_TRACK_PREFIX" "$sample_key" "$url" "$sample_key" "$assay_profile" "$factor" "$condition" >> "$trackdb"
done < "$SAMPLE_MANIFEST"

if is_true "$WRITE_IGV_SESSION"; then
    session="${OUTPUT_DIR}/09_browser/igv/cutnrun2tracks_session.xml"
    session_genome="$(awk -F '\t' 'NR==2 {print $5}' "$SAMPLE_MANIFEST")"
    printf '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<Session genome="%s" version="8">\n  <Resources>\n' "$session_genome" > "$session"
    find "${OUTPUT_DIR}/04_tracks" -type f -name '*.bw' -print | sort | while read -r bw; do
        printf '    <Resource path="%s"/>\n' "$bw" >> "$session"
    done
    printf '  </Resources>\n</Session>\n' >> "$session"
fi
run_logged python3 "${SCRIPT_DIR}/annotate_differential_results.py" "$OUTPUT_DIR"
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/07_annotation/stage_status.tsv"

#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

mkdir -p "${OUTPUT_DIR}/07_annotation" "${OUTPUT_DIR}/09_browser/ucsc" "${OUTPUT_DIR}/09_browser/igv"
trackdb="${OUTPUT_DIR}/09_browser/ucsc/trackDb.txt"
: > "$trackdb"

while IFS=$'\t' read -r cohort cohort_key genome assay_profile factor antibody_id layout target_class duplicate_policy caller peak_class n_samples sample_keys conditions; do
    [[ "$cohort" == "cohort_id" ]] && continue
    consensus_dir="${OUTPUT_DIR}/05_peaks/consensus/${cohort}/${caller}/${peak_class}"
    consensus="$(find "$consensus_dir" -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
    [[ -n "$consensus" ]] || continue
    out="${OUTPUT_DIR}/07_annotation/${cohort}/consensus"; mkdir -p "$out"
    if is_true "$RUN_SIMPLE_PEAK_ANNOTATION"; then
        gtf="$(reference_value GTF "$genome")"
        genes="${out}/genes.bed"
        awk 'BEGIN{FS=OFS="\t"} $3=="gene" {id="."; name="."; if(match($9,/gene_id "[^"]+"/))id=substr($9,RSTART+9,RLENGTH-10); if(match($9,/gene_name "[^"]+"/))name=substr($9,RSTART+11,RLENGTH-12); print $1,$4-1,$5,name,id,$7}' "$gtf" |
            sort -k1,1 -k2,2n > "$genes"
        bedtools closest -d -a "$consensus" -b "$genes" > "${out}/${cohort}.nearest_gene.tsv"
    fi
    if is_true "$RUN_CCRE_ANNOTATION"; then
        ccre="$(optional_reference_value CCRE_BED "$genome")"
        if [[ -n "$ccre" && "$ccre" != "." && -s "$ccre" ]]; then
            bedtools intersect -wao -a "$consensus" -b "$ccre" > "${out}/${cohort}.ccre_reference_overlaps.tsv"
        else
            printf '{"status":"SKIPPED","reason":"CCRE reference unavailable"}\n' > "${out}/ccre_SKIPPED.json"
        fi
    fi
done < "$COHORT_MANIFEST"

while IFS=$'\t' read -r sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition rest; do
    [[ "$sample_key" == "sample_key" ]] && continue
    bw="${OUTPUT_DIR}/04_tracks/cpm/${sample_key}.CPM.bw"
    [[ -s "$bw" ]] || continue
    url="$bw"
    [[ -n "$UCSC_BIGDATA_URL_BASE" ]] && url="${UCSC_BIGDATA_URL_BASE%/}/04_tracks/cpm/${sample_key}.CPM.bw"
    printf 'track %s_%s\ntype bigWig\nbigDataUrl %s\nshortLabel %s CPM\nlongLabel %s %s %s CPM\nvisibility full\nautoScale on\n\n' \
        "$UCSC_TRACK_PREFIX" "$sample_key" "$url" "$sample_key" "$ASSAY_PROFILE" "$factor" "$condition" >> "$trackdb"
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

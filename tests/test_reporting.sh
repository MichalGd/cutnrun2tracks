#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMPORARY="$(mktemp -d)"
trap 'rm -rf -- "$TEMPORARY"' EXIT
OUTPUT="${TEMPORARY}/output"
mkdir -p \
    "$OUTPUT/00_metadata" "$OUTPUT/03_alignment/metrics" \
    "$OUTPUT/04_tracks" "$OUTPUT/05_peaks/per_sample" "$OUTPUT/05_peaks/consensus" \
    "$OUTPUT/06_qc/alignment_and_complexity" "$OUTPUT/06_qc/frip_and_peak_reproducibility" \
    "$OUTPUT/08_differential" "$TEMPORARY/bin"

printf 'sample_key\tcohort_id\nTARGET.bioR1\tCOHORT_A\n' > "$OUTPUT/00_metadata/sample_manifest.tsv"
printf 'cohort_id\tn_samples\nCOHORT_A\t1\n' > "$OUTPUT/00_metadata/cohort_manifest.tsv"
printf 'sample_key\tlayout\tsignal_unit\tanalysis_observations\nTARGET.bioR1\tPE\tfragment\t1234\n' \
    > "$OUTPUT/06_qc/alignment_and_complexity/observation_counts.tsv"
printf 'sample_key\tlayout\tgenome\talignment_records\tspikein_mode\nTARGET.bioR1\tPE\thg38\t2000\tnone\n' \
    > "$OUTPUT/03_alignment/metrics/TARGET.bioR1.alignment.tsv"
printf 'sample_key\tq0_dup_retained\tq0_dup_removed\tq30_dup_retained\tq30_dup_removed\tanalysis_policy\nTARGET.bioR1\t1800\t1700\t1300\t1234\tremove\n' \
    > "$OUTPUT/03_alignment/metrics/TARGET.bioR1.filter_counts.tsv"
printf 'sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\tprimary_peak_count\tcaller_warnings\treason\nTARGET.bioR1\tCTRL.bioR1\tmacs3\tbroad\tSUCCESS\t42\t.\t.\n' \
    > "$OUTPUT/05_peaks/per_sample/peakcall_status.tsv"
printf 'cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\nCOHORT_A\tSUCCESS\t1\t1\t0\t42\t.\n' \
    > "$OUTPUT/05_peaks/consensus/consensus_status.tsv"
printf 'cohort_id\tpolicy\tstatus\treason\tlog\nCOHORT_A\tanalysis\tSUCCESS\t.\t.\n' \
    > "$OUTPUT/04_tracks/normalized_track_family_status.tsv"
printf 'status\tfailed_modules\tskipped_cohorts\nSUCCESS\t0\t0\n' > "$OUTPUT/08_differential/stage_status.tsv"
printf 'sample_key\tsignal_unit\ttotal\tin_consensus\tfrip\nTARGET.bioR1\tfragment\t1234\t321\t0.26012966\n' \
    > "$OUTPUT/06_qc/frip_and_peak_reproducibility/TARGET.bioR1.frip.tsv"

cat > "$TEMPORARY/bin/multiqc" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail
out=""; name=""
while (( $# )); do
    case "$1" in
        --outdir) out="$2"; shift 2 ;;
        --filename) name="$2"; shift 2 ;;
        --exclude|--ignore|--cl-config|--data-format) shift 2 ;;
        --export|--force) shift ;;
        *) shift ;;
    esac
done
[[ -n "$out" && -n "$name" ]]
mkdir -p "$out/${name}_data"
printf '<html><body>fake MultiQC report</body></html>\n' > "$out/${name}.html"
printf 'sample\tvalue\nTARGET.bioR1\t1\n' > "$out/${name}_data/multiqc_general_stats.txt"
echo 'MultiQC complete'
FAKE
chmod +x "$TEMPORARY/bin/multiqc"

before="$(sha256sum "$OUTPUT/00_metadata/sample_manifest.tsv")"
PATH="$TEMPORARY/bin:$PATH" bash "$ROOT/utilities/regenerate_reports.sh" --output-dir "$OUTPUT"
after="$(sha256sum "$OUTPUT/00_metadata/sample_manifest.tsv")"
[[ "$before" == "$after" ]]

for required in \
    pipeline_report.html run_summary.tsv warning_summary.tsv \
    cutnrun2tracks_multiqc_report.html multiqc_status.tsv \
    multiqc_custom_content_manifest.tsv report_checksums.sha256; do
    [[ -s "$OUTPUT/10_reports/$required" ]]
done
[[ -d "$OUTPUT/10_reports/cutnrun2tracks_multiqc_report_data" ]]
grep -q $'SUCCESS\t' "$OUTPUT/10_reports/multiqc_status.tsv"
grep -q 'cutnrun2tracks_observations' "$OUTPUT/10_reports/multiqc_custom_content_manifest.tsv"
grep -q 'Retained analysis observations' "$OUTPUT/10_reports/pipeline_report.html"

echo "Unified reporting and completed-run regeneration test passed"

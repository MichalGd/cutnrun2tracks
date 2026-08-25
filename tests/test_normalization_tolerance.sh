#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

mkdir -p "$TMP/bin" "$TMP/output/00_metadata"

cat > "$TMP/bin/Rscript" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
script="$1"; manifest="$2"; cohort="$3"; consensus="$4"; output="$5"; tables="$6"; policy="$7"
mkdir -p "$tables"
if [[ "$cohort" == "C_FAIL" ]]; then
    printf 'sample_key\tcohort_id\tpolicy\tconsensus_regions_with_any_count\tconsensus_count_sum\tstatus\n' \
        > "$tables/consensus_count_sums.tsv"
    printf 'FAIL_A.bioR1\tC_FAIL\tanalysis\t1\t0\tZERO\n' >> "$tables/consensus_count_sums.tsv"
    printf 'FAIL_B.bioR2\tC_FAIL\tanalysis\t1\t4\tNONZERO\n' >> "$tables/consensus_count_sums.tsv"
    echo 'zero consensus counts for samples: FAIL_A.bioR1' >&2
    exit 42
fi
printf 'sample_key\tcohort_id\tpolicy\tsignal_unit\tsize_factor\tdeseq2_consensus_scale\tconsensus_count_sum\tcohort_geometric_mean_column_sum\trobust_effective_library_size\tdeseq2_robust_cpm_scale\n' \
    > "$tables/normalization_factors.tsv"
printf 'OK_A.bioR1\tC_OK\tanalysis\tfragment\t1\t1\t10\t10\t10\t100000\n' >> "$tables/normalization_factors.tsv"
printf 'OK_B.bioR2\tC_OK\tanalysis\tfragment\t1\t1\t10\t10\t10\t100000\n' >> "$tables/normalization_factors.tsv"
printf 'region_id\tOK_A.bioR1\tOK_B.bioR2\nregion1\t10\t10\n' > "$tables/raw_counts.tsv.gz"
MOCK
chmod +x "$TMP/bin/Rscript"

header='sample_key\tis_control\tcohort_id'
printf '%b\n' "$header" > "$TMP/output/00_metadata/sample_manifest.tsv"
printf 'FAIL_A.bioR1\tFALSE\tC_FAIL\nFAIL_B.bioR2\tFALSE\tC_FAIL\nOK_A.bioR1\tFALSE\tC_OK\nOK_B.bioR2\tFALSE\tC_OK\n' \
    >> "$TMP/output/00_metadata/sample_manifest.tsv"

cohort_header='cohort_id\tcohort_key\tgenome\tassay_profile\tfactor\tantibody_id\tlayout\ttarget_class\tanalysis_duplicate_policy\tprimary_peak_caller\tprimary_peak_class\tn_biological_samples\tsample_keys\tconditions'
printf '%b\n' "$cohort_header" > "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'C_FAIL\tkey1\thg38\tcuttag\tF1\tA1\tPE\tbroad\tretain\tmacs3\tbroad\t2\tFAIL_A.bioR1,FAIL_B.bioR2\tA,B\n' \
    >> "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'C_OK\tkey2\thg38\tcuttag\tF2\tA2\tPE\tbroad\tretain\tmacs3\tbroad\t2\tOK_A.bioR1,OK_B.bioR2\tA,B\n' \
    >> "$TMP/output/00_metadata/cohort_manifest.tsv"

for cohort in C_FAIL C_OK; do
    directory="$TMP/output/05_peaks/consensus/$cohort/macs3/broad"
    mkdir -p "$directory"
    printf 'chr1\t0\t100\t%s.region1\t2\n' "$cohort" > "$directory/$cohort.consensus.bed"
done

cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR=$TMP/output
REQUIRE_ALL_ENABLED_TRACKS=false
GENERATE_DESEQ2_CONSENSUS_TRACKS=false
RUN_DIFFBIND=false
RUN_DESEQ2_ENRICHMENT=true
GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS=false
GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS=false
GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS=false
RUN_CONTROL_SUBTRACTED_SENSITIVITY=false
RUN_TARGET_CONTROL_INTERACTION=false
DIFFERENTIAL_NORMALIZATION=internal
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/normalized_tracks_batch.sh"
grep -q '^COMPLETED_WITH_WARNINGS' "$TMP/output/04_tracks/stage_status.tsv"
grep -q 'zero consensus counts for samples: FAIL_A.bioR1' \
    "$TMP/output/logs/normalized_tracks/C_FAIL.analysis.factors.log"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_FAIL/SKIPPED.json"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_FAIL/tables/consensus_count_sums.tsv"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_OK/tables/raw_counts.tsv.gz"
awk -F '\t' '$1=="C_FAIL" {found=($3=="SKIPPED")} END {exit !found}' \
    "$TMP/output/04_tracks/normalized_track_family_status.tsv"
awk -F '\t' '$1=="C_OK" {found=($3=="SUCCESS")} END {exit !found}' \
    "$TMP/output/04_tracks/normalized_track_family_status.tsv"

sed 's/RUN_DESEQ2_ENRICHMENT=true/RUN_DESEQ2_ENRICHMENT=false/' "$TMP/config.sh" > "$TMP/config.differential.sh"
C2T_CONFIG="$TMP/config.differential.sh" bash "$ROOT/scripts/differential_batch.sh"
grep -q '^COMPLETED_WITH_WARNINGS' "$TMP/output/08_differential/stage_status.tsv"
test -s "$TMP/output/08_differential/C_FAIL/broad/SKIPPED.json"
test ! -e "$TMP/output/08_differential/C_FAIL/broad/FAILED.json"

sed 's/REQUIRE_ALL_ENABLED_TRACKS=false/REQUIRE_ALL_ENABLED_TRACKS=true/' "$TMP/config.sh" > "$TMP/config.strict.sh"
if PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.strict.sh" \
    bash "$ROOT/scripts/normalized_tracks_batch.sh" >/dev/null 2>&1; then
    echo "ERROR: strict normalized-track policy accepted a zero-count cohort" >&2
    exit 1
fi
grep -q '^FAILED' "$TMP/output/04_tracks/stage_status.tsv"

echo "Normalized-track cohort continuation and strict-policy regression tests passed"

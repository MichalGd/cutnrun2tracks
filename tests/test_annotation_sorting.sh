#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

mkdir -p "$TMP/bin" "$TMP/output/00_metadata" \
    "$TMP/output/05_peaks/consensus/C1/macs3/broad"

cat > "$TMP/bin/bedtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
command="$1"; shift
case "$command" in
    sort)
        genome=""; input=""
        while (($#)); do
            case "$1" in
                -faidx) genome="$2"; shift 2 ;;
                -i) input="$2"; shift 2 ;;
                *) exit 2 ;;
            esac
        done
        printf 'sort\t%s\t%s\n' "$genome" "$input" >> "$MOCK_BEDTOOLS_LOG"
        awk 'NR==FNR {rank[$1]=NR; next} {print rank[$1],$0}' OFS='\t' "$genome" "$input" |
            sort -k1,1n -k3,3n | cut -f2-
        ;;
    closest)
        genome=""; a=""; b=""
        while (($#)); do
            case "$1" in
                -a) a="$2"; shift 2 ;;
                -b) b="$2"; shift 2 ;;
                -g) genome="$2"; shift 2 ;;
                -d) shift ;;
                *) exit 2 ;;
            esac
        done
        [[ -n "$genome" && -s "$genome" && -s "$a" && -s "$b" ]]
        printf 'closest\t%s\t%s\t%s\n' "$genome" "$a" "$b" >> "$MOCK_BEDTOOLS_LOG"
        cat "$a"
        ;;
    *) exit 2 ;;
esac
MOCK

cat > "$TMP/bin/python3" <<'MOCK'
#!/usr/bin/env bash
exit 0
MOCK
chmod +x "$TMP/bin/bedtools" "$TMP/bin/python3"

printf 'chr1\t1000\nchr2\t1000\n' > "$TMP/chrom.sizes"
printf 'chr2\tsource\tgene\t21\t40\t.\t+\t.\tgene_id "G2"; gene_name "Gene2";\n' > "$TMP/annotation.gtf"
printf 'chr1\tsource\tgene\t11\t30\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";\n' >> "$TMP/annotation.gtf"
printf 'chr2\t20\t40\tC1.region2\t2\nchr1\t10\t30\tC1.region1\t2\n' \
    > "$TMP/output/05_peaks/consensus/C1/macs3/broad/C1.consensus.bed"

printf 'cohort_id\tcohort_key\tgenome\tassay_profile\tfactor\tantibody_id\tlayout\ttarget_class\tanalysis_duplicate_policy\tprimary_peak_caller\tprimary_peak_class\tn_biological_samples\tsample_keys\tconditions\n' \
    > "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'C1\tkey\thg38\tcuttag\tH3K27ac\tAB1\tPE\tbroad\tretain\tmacs3\tbroad\t2\tS1,S2\tA,B\n' \
    >> "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'sample_key\tsample_id\treplicate\tlayout\tgenome\tassay_profile\tfactor\tantibody_id\ttarget_class\tcondition\trest\n' \
    > "$TMP/output/00_metadata/sample_manifest.tsv"

cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR=$TMP/output
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=false
WRITE_IGV_SESSION=false
GTF_HG38=$TMP/annotation.gtf
CHROM_SIZES_HG38=$TMP/chrom.sizes
UCSC_BIGDATA_URL_BASE=
UCSC_TRACK_PREFIX=test
ASSAY_PROFILE=cuttag
EOF

MOCK_BEDTOOLS_LOG="$TMP/bedtools.log" PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" \
    bash "$ROOT/scripts/annotate_browser.sh"

test "$(grep -c '^sort' "$TMP/bedtools.log")" -eq 2
grep -q $'^closest\t.*/chrom.sizes\t' "$TMP/bedtools.log"
test "$(awk 'NR==1 {print $1}' "$TMP/output/07_annotation/C1/consensus/genes.bed")" = chr1
test "$(awk 'NR==1 {print $1}' "$TMP/output/07_annotation/C1/consensus/C1.consensus.sorted.bed")" = chr1
grep -qx 'SUCCESS' "$TMP/output/07_annotation/stage_status.tsv"

echo "Genome-order-safe annotation regression test passed"

#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT
mkdir -p "$TMP/bin" "$TMP/output/00_metadata" "$TMP/output/03_alignment/analysis" "$TMP/output/logs"

cat > "$TMP/bin/samtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
    view)
        shift
        output=""
        while (($#)); do
            if [[ "$1" == "-o" ]]; then output="$2"; shift 2; else shift; fi
        done
        if [[ -n "$output" ]]; then cat > "$output"; else
            printf '@HD\tVN:1.6\nread1\t99\tchr1\t1\t60\t50M\t=\t51\t100\tA\tI\n'
        fi
        ;;
    sort)
        shift
        output="" input=""
        while (($#)); do
            case "$1" in -@) shift 2 ;; -o) output="$2"; shift 2 ;; *) input="$1"; shift ;; esac
        done
        cp "$input" "$output"
        ;;
    *) exit 2 ;;
esac
MOCK

cat > "$TMP/bin/bedtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
[[ "$1" == "genomecov" ]] || exit 2
printf 'chr1\t0\t100\t1\n'
MOCK

cat > "$TMP/bin/mock_macs3" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
out="" name=""
while (($#)); do
    case "$1" in --outdir) out="$2"; shift 2 ;; -n) name="$2"; shift 2 ;; *) shift ;; esac
done
mkdir -p "$out"
: > "$out/${name}_peaks.broadPeak"
MOCK

cat > "$TMP/bin/mock_seacr" <<'MOCK'
#!/usr/bin/env bash
exit 42
MOCK
chmod +x "$TMP/bin/"*

header='sample_key\tsample_id\treplicate\tlayout\tgenome\tassay_profile\tfactor\tantibody_id\ttarget_class\tcondition\ttreatment\tcell_type\tis_control\tcontrol_type\tcontrol_id\tcontrol_key\tanalysis_duplicate_policy\tblacklist\tspikein_to_host_ratio\tspikein_stage\tspikein_lot\tbatch\tdonor\toutput_prefix\ttechnical_units\tfastq_1_list\tfastq_2_list\tcohort_id\tcohort_key\tprimary_peak_caller\tprimary_peak_class'
target='TARGET.bioR1\tTARGET\t1\tPE\thg38\tcuttag\tH3K27ac\tAB1\tbroad\tA\tnone\tcell\tFALSE\tnone\tCTRL\tCTRL.bioR1\tretain\t.\t.\t.\t.\tb1\td1\tTARGET\t1\ta\tb\tC1\tkey\tmacs3\tbroad'
control='CTRL.bioR1\tCTRL\t1\tPE\thg38\tcuttag\tIgG\tIgG\tnarrow\tA\tnone\tcell\tTRUE\tigg\tCTRL\t.\tretain\t.\t.\t.\t.\tb1\td1\tCTRL\t1\tc\td\t.\t.\tnone\tcontrol'
printf '%b\n%b\n%b\n' "$header" "$target" "$control" > "$TMP/output/00_metadata/sample_manifest.tsv"
printf 'chr1\t1000\n' > "$TMP/chrom.sizes"
touch "$TMP/output/03_alignment/analysis/TARGET.bioR1.host.analysis.bam"
touch "$TMP/output/03_alignment/analysis/CTRL.bioR1.host.analysis.bam"

cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR=$TMP/output
PEAK_CALLERS=seacr,macs3
PEAKCALL_FAILURE_POLICY=continue
PEAKCALL_PARALLEL_JOBS=1
ALLOW_EMPTY_PEAKS=true
MACS3_COMMAND=$TMP/bin/mock_macs3
MACS3_GENOME_SIZE=auto
MACS3_QVALUE=0.01
MACS3_KEEP_DUP=all
MACS3_CALL_SUMMITS=false
MACS3_BROAD_CUTOFF=0.1
MACS3_CUTRUN_SE_SHIFT=-75
MACS3_CUTRUN_SE_EXTSIZE=150
MACS3_CUTTAG_SE_SHIFT=-75
MACS3_CUTTAG_SE_EXTSIZE=150
SEACR_COMMAND=$TMP/bin/mock_seacr
SEACR_MAX_FRAGMENT=1000
SEACR_NO_CONTROL_THRESHOLD=0.05
SEACR_CONTROL_NORMALIZATION=norm
SEACR_MODE=stringent
THREADS_SAMTOOLS=1
CHROM_SIZES_HG38=$TMP/chrom.sizes
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/peakcall_batch.sh"
grep -q '^COMPLETED_WITH_WARNINGS' "$TMP/output/05_peaks/per_sample/stage_status.tsv"
awk -F '\t' 'NR==2 {exit !($5=="EMPTY" && $7 ~ /macs3:broad=EMPTY/ && $7 ~ /seacr:narrow=ERROR/)}' \
    "$TMP/output/05_peaks/per_sample/TARGET.bioR1/peakcall_metadata.tsv"

sed 's/PEAKCALL_FAILURE_POLICY=continue/PEAKCALL_FAILURE_POLICY=fail/' "$TMP/config.sh" > "$TMP/config.strict.sh"
if PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.strict.sh" bash "$ROOT/scripts/peakcall_batch.sh" >/dev/null 2>&1; then
    echo "ERROR: strict peak-call policy accepted failed callers" >&2
    exit 1
fi
echo "Peak-caller continuation and strict-policy regression tests passed"

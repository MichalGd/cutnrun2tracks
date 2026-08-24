#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temporary="$(mktemp -d)"
trap 'rm -rf -- "$temporary"' EXIT

fake_bin="$temporary/bin"
output="$temporary/output"
mkdir -p "$fake_bin" "$output/00_metadata" \
    "$output/03_alignment/sorted" "$output/03_alignment/marked" \
    "$output/03_alignment/metrics"

cat > "$fake_bin/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
command_name="${1:?command required}"
shift
case "$command_name" in
    quickcheck|idxstats)
        exit 0
        ;;
    view)
        count=false
        destination=""
        while (( $# > 0 )); do
            case "$1" in
                -c) count=true; shift ;;
                -o) destination="$2"; shift 2 ;;
                -@|-q|-f|-F) shift 2 ;;
                -b) shift ;;
                *) shift ;;
            esac
        done
        if [[ "$count" == "true" ]]; then
            printf '1\n'
        else
            [[ -n "$destination" ]]
            printf 'synthetic bam\n' > "$destination"
        fi
        ;;
    sort)
        destination=""
        while (( $# > 0 )); do
            case "$1" in
                -o) destination="$2"; shift 2 ;;
                -@) shift 2 ;;
                -n) shift ;;
                *) shift ;;
            esac
        done
        [[ -n "$destination" ]]
        printf 'synthetic sorted bam\n' > "$destination"
        ;;
    fixmate)
        destination="${*: -1}"
        printf 'synthetic fixmate bam\n' > "$destination"
        ;;
    index)
        bam="${*: -1}"
        printf 'synthetic index\n' > "${bam}.bai"
        ;;
    *)
        echo "unexpected fake samtools command: $command_name" >&2
        exit 97
        ;;
esac
EOF

cat > "$fake_bin/bedtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" == "intersect" ]]
shift
input=""
while (( $# > 0 )); do
    case "$1" in
        -abam) input="$2"; shift 2 ;;
        -b) shift 2 ;;
        *) shift ;;
    esac
done
[[ -n "$input" ]]
cat "$input"
EOF

cat > "$fake_bin/picard" <<'EOF'
#!/usr/bin/env bash
echo "ERROR: Picard should not run when the marked BAM is reusable" >&2
exit 98
EOF

chmod +x "$fake_bin/samtools" "$fake_bin/bedtools" "$fake_bin/picard"

key="sample.bioR1"
printf 'synthetic sorted bam\n' > "$output/03_alignment/sorted/${key}.host.sorted.bam"
printf 'synthetic marked bam\n' > "$output/03_alignment/marked/${key}.host.marked.bam"
printf 'synthetic marked index\n' > "$output/03_alignment/marked/${key}.host.marked.bam.bai"
printf 'synthetic metrics\n' > "$output/03_alignment/metrics/${key}.duplicate_metrics.txt"

canonical="$temporary/hg38.canonical_contigs.txt"
blacklist="$temporary/hg38.blacklist.bed"
printf 'chr1\n' > "$canonical"
: > "$blacklist"

manifest="$output/00_metadata/sample_manifest.tsv"
printf '%s\n' \
    $'sample_key\tsample_id\treplicate\tlayout\tgenome\tassay_profile\tfactor\tantibody_id\ttarget_class\tcondition\ttreatment\tcell_type\tis_control\tcontrol_type\tcontrol_id\tcontrol_key\tanalysis_duplicate_policy\tblacklist\trest' \
    "$key"$'\tsample\t1\tPE\thg38\tcuttag\tH3K27ac\tH3K27ac\tbroad\ttreated\tnone\tcells\tFALSE\tnone\tcontrol\tcontrol.bioR1\tretain\t'"$blacklist"$'\t.' \
    > "$manifest"

config="$temporary/config.conf"
cat > "$config" <<EOF
OUTPUT_DIR=$output
THREADS_PARALLEL_JOBS=1
THREADS_SAMTOOLS=1
REMOVE_MITO=false
ALLOW_EMPTY_FILTERED_BAM=false
PERMISSIVE_MIN_MAPQ=0
INTERMEDIATE_MIN_MAPQ=0
MIN_MAPQ=30
PICARD_COMMAND=picard
CANONICAL_CONTIGS_HG38=$canonical
EOF

log="$temporary/filtering.log"
PATH="$fake_bin:$PATH" C2T_CONFIG="$config" \
    bash "$ROOT/scripts/mark_filter_batch.sh" > "$log" 2>&1

grep -qx 'SUCCESS' "$output/03_alignment/analysis/stage_status.tsv"
grep -q "Reusing validated marked BAM for $key" \
    "$output/logs/filtering/${key}.picard.log"

for branch in q0_dup-retained q0_dup-removed q30_dup-retained q30_dup-removed; do
    count="$(find "$output/03_alignment/filtered/$branch" -maxdepth 1 -type f -name '*.bam' | wc -l)"
    [[ "$count" -eq 1 ]]
done

if find "$output/03_alignment/filtered" -maxdepth 1 -type d -name '.filter.*' -print -quit | grep -q .; then
    echo "ERROR: filtering temporary directory was not cleaned" >&2
    exit 1
fi
if grep -q 'tmp: unbound variable' "$log"; then
    echo "ERROR: filtering cleanup trap leaked into the worker" >&2
    exit 1
fi

echo "Filtering cleanup-scope regression test passed"

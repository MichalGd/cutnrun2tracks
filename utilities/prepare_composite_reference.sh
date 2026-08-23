#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: prepare_composite_reference.sh <host.fa> <spike.fa> <spike_reference_id> <output_dir> [threads]" >&2
}
(( $# >= 4 && $# <= 5 )) || { usage; exit 2; }
HOST_FASTA="$1"; SPIKE_FASTA="$2"; REFERENCE_ID="$3"; OUTPUT_DIR="$4"; THREADS="${5:-4}"
[[ -s "$HOST_FASTA" && -s "$SPIKE_FASTA" ]] || { echo "ERROR: input FASTA missing" >&2; exit 1; }
[[ "$REFERENCE_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "ERROR: unsafe spike reference ID" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR"
resolved="$(realpath -m "$OUTPUT_DIR")"
[[ "$resolved" != "/" && ${#resolved} -gt 8 ]] || { echo "ERROR: unsafe output directory" >&2; exit 1; }
spike_namespaced="$OUTPUT_DIR/${REFERENCE_ID}.namespaced.fa"
composite="$OUTPUT_DIR/host_plus_${REFERENCE_ID}.fa"
allowed="$OUTPUT_DIR/${REFERENCE_ID}.allowed_contigs.txt"

awk -v prefix="${REFERENCE_ID}__" '/^>/ {$0=">"prefix substr($0,2)} {print}' "$SPIKE_FASTA" > "$spike_namespaced"
grep '^>' "$spike_namespaced" | sed 's/^>//;s/[[:space:]].*$//' > "$allowed"
host_headers="$(mktemp)"; spike_headers="$(mktemp)"; trap 'rm -f "$host_headers" "$spike_headers"' EXIT
grep '^>' "$HOST_FASTA" | sed 's/^>//;s/[[:space:]].*$//' | sort -u > "$host_headers"
sort -u "$allowed" > "$spike_headers"
if comm -12 "$host_headers" "$spike_headers" | grep -q .; then
    echo "ERROR: contig collision after namespacing" >&2
    exit 1
fi
cp "$HOST_FASTA" "$composite"
dd if="$spike_namespaced" of="$composite" oflag=append conv=notrunc status=none
samtools faidx "$composite"
bowtie2-build --threads "$THREADS" "$composite" "$OUTPUT_DIR/host_plus_${REFERENCE_ID}"
{
    printf 'item\tpath\tsha256\n'
    for file in "$HOST_FASTA" "$SPIKE_FASTA" "$spike_namespaced" "$composite" "$allowed"; do
        printf '%s\t%s\t%s\n' "$(basename "$file")" "$(realpath "$file")" "$(sha256sum "$file" | awk '{print $1}')"
    done
} > "$OUTPUT_DIR/reference_manifest.tsv"
echo "Composite reference created: $OUTPUT_DIR/host_plus_${REFERENCE_ID}"

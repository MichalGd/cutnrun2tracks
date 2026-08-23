#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/03_alignment/sorted" "${OUTPUT_DIR}/03_alignment/spikein/composite" \
    "${OUTPUT_DIR}/03_alignment/spikein/spike" "${OUTPUT_DIR}/03_alignment/metrics" "${OUTPUT_DIR}/logs/alignment"

worker() {
    local key="$1" layout="$2" genome="$3"
    local fq1="${OUTPUT_DIR}/02_trimmed_fastq/${key}.R1.trimmed.fastq.gz"
    local fq2="${OUTPUT_DIR}/02_trimmed_fastq/${key}.R2.trimmed.fastq.gz"
    local host="${OUTPUT_DIR}/03_alignment/sorted/${key}.host.sorted.bam"
    local composite="${OUTPUT_DIR}/03_alignment/spikein/composite/${key}.composite.sorted.bam"
    local target="$host" index log="${OUTPUT_DIR}/logs/alignment/${key}.bowtie2.log"
    [[ -s "$fq1" ]] || die "trimmed R1 missing for $key"
    if [[ "$SPIKEIN_MODE" == "none" ]]; then
        index="$(reference_value INDEX "$genome")"
    else
        index="$SPIKEIN_INDEX"
        target="$composite"
    fi
    local args=(--"$BOWTIE2_MODE" --"$BOWTIE2_PRESET" --seed "$BOWTIE2_SEED" --rg-id "$key" --rg "SM:$key")
    [[ "$BOWTIE2_DOVETAIL" == "true" && "$layout" == "PE" ]] && args+=(--dovetail)
    [[ "$BOWTIE2_MIXED" == "false" && "$layout" == "PE" ]] && args+=(--no-mixed)
    [[ "$BOWTIE2_DISCORDANT" == "false" && "$layout" == "PE" ]] && args+=(--no-discordant)
    [[ "$layout" == "PE" ]] && args+=(-I "$BOWTIE2_MIN_INSERT" -X "$BOWTIE2_MAX_INSERT")
    local extra=()
    [[ -n "$BOWTIE2_EXTRA_ARGS" ]] && read -r -a extra <<< "$BOWTIE2_EXTRA_ARGS"
    record_command bowtie2 "${args[@]}" "${extra[@]}" -x "$index"
    if [[ "$layout" == "PE" ]]; then
        bowtie2 "${args[@]}" "${extra[@]}" -p "$THREADS_BOWTIE2" -x "$index" -1 "$fq1" -2 "$fq2" 2>"$log" |
            samtools sort -@ "$THREADS_SAMTOOLS" -o "$target" -
    else
        bowtie2 "${args[@]}" "${extra[@]}" -p "$THREADS_BOWTIE2" -x "$index" -U "$fq1" 2>"$log" |
            samtools sort -@ "$THREADS_SAMTOOLS" -o "$target" -
    fi
    samtools quickcheck "$target"
    samtools index -@ "$THREADS_SAMTOOLS" "$target"
    if [[ "$SPIKEIN_MODE" != "none" ]]; then
        local chrom_sizes host_contigs=() spike_contigs=()
        chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
        mapfile -t host_contigs < <(awk 'NF>=2 {print $1}' "$chrom_sizes")
        mapfile -t spike_contigs < <(awk 'NF {print $1}' "$SPIKEIN_ALLOWED_CONTIGS")
        samtools view -@ "$THREADS_SAMTOOLS" -b -o "$host" "$composite" "${host_contigs[@]}"
        samtools view -@ "$THREADS_SAMTOOLS" -b -o "${OUTPUT_DIR}/03_alignment/spikein/spike/${key}.${SPIKEIN_MODE}.bam" \
            "$composite" "${spike_contigs[@]}"
        samtools index -@ "$THREADS_SAMTOOLS" "$host"
        samtools index -@ "$THREADS_SAMTOOLS" "${OUTPUT_DIR}/03_alignment/spikein/spike/${key}.${SPIKEIN_MODE}.bam"
    fi
    local records
    records="$(samtools view -c "$host")"
    printf 'sample_key\tlayout\tgenome\talignment_records\tspikein_mode\n%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$layout" "$genome" "$records" "$SPIKEIN_MODE" > "${OUTPUT_DIR}/03_alignment/metrics/${key}.alignment.tsv"
}

parallel_pool_init "$THREADS_PARALLEL_JOBS"
while IFS=$'\t' read -r sample_key sample_id replicate layout genome rest; do
    [[ "$sample_key" == "sample_key" ]] && continue
    parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all
printf 'status\nSUCCESS\n' > "${OUTPUT_DIR}/03_alignment/sorted/stage_status.tsv"

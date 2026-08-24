#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

required=(python3 bowtie2 samtools bedtools "${PICARD_COMMAND}" "${MACS3_COMMAND}" bamCoverage Rscript)
is_true "$TRIM_ADAPTERS" && required+=(trim_galore)
is_true "$RUN_FASTQC" && required+=(fastqc)
is_true "$RUN_MULTIQC" && required+=(multiqc)
is_true "$GENERATE_COVERAGE_BIGWIGS" && required+=(bedGraphToBigWig)
[[ ",$PEAK_CALLERS," == *,seacr,* ]] && required+=("${SEACR_COMMAND}")
is_true "$RUN_PRESEQ" && required+=(preseq)
is_true "$RUN_ATAQV_QC" && required+=(ataqv)
is_true "$GENERATE_ATAQV_VIEWER" && required+=(mkarv)
is_true "$RUN_TSS_SIGNAL_PROFILE" && required+=(computeMatrix plotProfile)
is_true "$RUN_METAGENE" && required+=(computeMatrix plotProfile plotHeatmap)
required+=(plotFingerprint)
is_true "$RUN_MOTIF_ENRICHMENT" && die "RUN_MOTIF_ENRICHMENT is not implemented in this release"

missing=()
for command_name in "${required[@]}"; do
    if [[ "$command_name" == */* ]]; then
        [[ -x "$command_name" || -f "$command_name" ]] || missing+=("$command_name")
    elif ! command -v "$command_name" >/dev/null 2>&1; then
        missing+=("$command_name")
    fi
done
(( ${#missing[@]} == 0 )) || die "missing required commands: ${missing[*]}"

if is_true "$RUN_METAGENE"; then
    [[ -s "$METAGENE_GENE_SET_MANIFEST" ]] || die "metagene gene-set manifest missing: $METAGENE_GENE_SET_MANIFEST"
    python3 -c 'import pyBigWig' || die "RUN_METAGENE requires the pyBigWig Python package"
fi

mkdir -p "${OUTPUT_DIR}/00_metadata"
versions="${OUTPUT_DIR}/00_metadata/software_versions.tsv"
printf 'tool\tversion\n' > "$versions"
version_commands=(python3 bowtie2 samtools bedtools bamCoverage Rscript)
is_true "$RUN_METAGENE" && version_commands+=(computeMatrix plotProfile plotHeatmap)
for command_name in "${version_commands[@]}"; do
    version="$($command_name --version 2>&1 | head -n 1 || true)"
    printf '%s\t%s\n' "$command_name" "${version//$'\t'/ }" >> "$versions"
done

tail -n +2 "$SAMPLE_MANIFEST" | while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    index="$(reference_value INDEX "$genome")"
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    canonical_contigs="$(reference_value CANONICAL_CONTIGS "$genome")"
    fasta="$(reference_value FASTA "$genome")"
    gtf="$(reference_value GTF "$genome")"
    [[ -s "$chrom_sizes" && -s "$canonical_contigs" && -s "$fasta" && -s "$gtf" && -s "$blacklist" ]] || \
        die "missing reference file for $sample_key"
    if is_true "$RUN_CCRE_ANNOTATION"; then
        ccre="$(optional_reference_value CCRE_BED "$genome")"
        [[ -s "$ccre" ]] || die "RUN_CCRE_ANNOTATION requires CCRE_BED for $genome"
    fi
    [[ -f "${index}.1.bt2" || -f "${index}.1.bt2l" ]] || die "Bowtie2 index prefix invalid: $index"
done

if [[ "$SPIKEIN_MODE" != "none" ]]; then
    [[ -s "$SPIKEIN_FASTA" && -s "$SPIKEIN_CHROM_SIZES" && -s "$SPIKEIN_ALLOWED_CONTIGS" ]] || \
        die "spike-in reference manifest files are incomplete"
    [[ -f "${SPIKEIN_INDEX}.1.bt2" || -f "${SPIKEIN_INDEX}.1.bt2l" ]] || \
        die "composite Bowtie2 index prefix invalid: $SPIKEIN_INDEX"
fi

printf 'status\tprofile\tspikein_mode\tpeak_callers\nSUCCESS\t%s\t%s\t%s\n' \
    "$ASSAY_PROFILE" "$SPIKEIN_MODE" "$PEAK_CALLERS" > "${OUTPUT_DIR}/00_metadata/preflight_status.tsv"
note "Preflight successful"

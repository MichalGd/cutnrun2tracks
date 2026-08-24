#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

qc_root="${OUTPUT_DIR}/06_qc"
mkdir -p "$qc_root/alignment_and_complexity" "$qc_root/fragment_length_and_periodicity" \
    "$qc_root/frip_and_peak_reproducibility" "$qc_root/correlation_pca_fingerprint" \
    "$qc_root/tss_signal_profile" "$qc_root/controls" "$qc_root/experimental_ATAC_derived_ataqv" \
    "${OUTPUT_DIR}/logs/qc"
summary="$qc_root/alignment_and_complexity/observation_counts.tsv"
printf 'sample_key\tlayout\tsignal_unit\tanalysis_observations\n' > "$summary"

while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" ]] && continue
    bam="$(analysis_bam_path "$sample_key")"
    count="$(signal_count "$bam" "$layout")"
    unit="$([[ "$layout" == "PE" ]] && echo fragment || echo read)"
    printf '%s\t%s\t%s\t%s\n' "$sample_key" "$layout" "$unit" "$count" >> "$summary"
    samtools flagstat -@ "$THREADS_SAMTOOLS" "$bam" > "$qc_root/alignment_and_complexity/${sample_key}.flagstat.txt"
    samtools stats -@ "$THREADS_SAMTOOLS" "$bam" > "$qc_root/alignment_and_complexity/${sample_key}.samtools_stats.txt"
    if [[ "$layout" == "PE" ]] && is_true "$RUN_FRAGMENT_QC"; then
        samtools view -f 66 -F 3840 "$bam" |
            awk -v maximum="$FRAGMENT_PLOT_MAX_BP" 'BEGIN{OFS="\t"} {t=$9; if(t<0)t=-t; if(t>0 && t<=maximum)n[t]++} END{print "fragment_length","count"; for(i=1;i<=maximum;i++)print i,n[i]+0}' \
            > "$qc_root/fragment_length_and_periodicity/${sample_key}.fragment_lengths.tsv"
    fi
    if is_true "$RUN_PRESEQ"; then
        preseq lc_extrap -B -o "$qc_root/alignment_and_complexity/${sample_key}.preseq.txt" "$bam" \
            >"${OUTPUT_DIR}/logs/qc/${sample_key}.preseq.log" 2>&1 || warn "preseq failed for $sample_key"
    fi
    if [[ "$is_control" == "FALSE" ]]; then
        consensus_dir="${OUTPUT_DIR}/05_peaks/consensus/${cohort_id}/${primary_caller}/${primary_class}"
        consensus="$(find "$consensus_dir" -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
        if [[ -n "$consensus" ]]; then
            if [[ "$layout" == "PE" ]]; then
                tmp="$(mktemp -d "$qc_root/frip_and_peak_reproducibility/.frip.XXXXXX")"
                samtools sort -n -@ "$THREADS_SAMTOOLS" -o "$tmp/name.bam" "$bam"
                bedtools bamtobed -bedpe -i "$tmp/name.bam" |
                    awk 'BEGIN{OFS="\t"} $1==$4 {start=($2<$5?$2:$5); end=($3>$6?$3:$6); if(end>start)print $1,start,end}' > "$tmp/fragments.bed"
                total="$(wc -l < "$tmp/fragments.bed")"
                in_peaks="$(bedtools intersect -u -a "$tmp/fragments.bed" -b "$consensus" | wc -l)"
                rm -rf "$tmp"
            else
                total="$count"
                in_peaks="$(bedtools intersect -u -abam "$bam" -b "$consensus" | samtools view -c -)"
            fi
            frip="$(awk -v a="$in_peaks" -v n="$total" 'BEGIN{if(n>0)printf "%.8f",a/n; else print "NA"}')"
            printf 'sample_key\tsignal_unit\ttotal\tin_consensus\tfrip\n%s\t%s\t%s\t%s\t%s\n' \
                "$sample_key" "$unit" "$total" "$in_peaks" "$frip" \
                > "$qc_root/frip_and_peak_reproducibility/${sample_key}.frip.tsv"
        fi
        if [[ "$control_key" != "." ]]; then
            control_bam="$(analysis_bam_path "$control_key")"
            fingerprint_args=(-b "$bam" "$control_bam" --labels "$sample_key" "$control_key"
                --plotFile "$qc_root/controls/${sample_key}.target_control_fingerprint.png"
                --outRawCounts "$qc_root/controls/${sample_key}.target_control_fingerprint.tsv"
                --numberOfProcessors "$THREADS_BAMCOVERAGE")
            [[ "$layout" == "PE" ]] && fingerprint_args+=(--samFlagInclude 66 --extendReads)
            plotFingerprint "${fingerprint_args[@]}" >"${OUTPUT_DIR}/logs/qc/${sample_key}.fingerprint.log" 2>&1 || \
                warn "fingerprint failed for $sample_key"
        fi
    fi
done < "$SAMPLE_MANIFEST"

if is_true "$RUN_TSS_SIGNAL_PROFILE"; then
    while IFS=$'\t' read -r sample_key sample_id replicate layout genome rest; do
        [[ "$sample_key" == "sample_key" ]] && continue
        tss="$(optional_reference_value TSS_BED "$genome")"
        if [[ -z "$tss" || "$tss" == "." ]]; then
            tss="$qc_root/tss_signal_profile/reference/${genome}.tss.bed"
            if [[ ! -s "$tss" ]]; then
                gtf="$(reference_value GTF "$genome")"
                python3 "${SCRIPT_DIR}/prepare_tss_bed.py" "$gtf" "$tss"
            fi
        fi
        bw="${OUTPUT_DIR}/04_tracks/cpm/${sample_key}.CPM.bw"
        [[ -s "$bw" ]] || continue
        computeMatrix reference-point --referencePoint TSS -b "$TSS_PROFILE_UPSTREAM" -a "$TSS_PROFILE_DOWNSTREAM" \
            -R "$tss" -S "$bw" -o "$qc_root/tss_signal_profile/${sample_key}.matrix.gz" \
            --numberOfProcessors "$THREADS_BAMCOVERAGE"
        plotProfile -m "$qc_root/tss_signal_profile/${sample_key}.matrix.gz" \
            -out "$qc_root/tss_signal_profile/${sample_key}.descriptive_TSS_profile.png" --plotTitle "$sample_key descriptive TSS signal"
    done < "$SAMPLE_MANIFEST"
fi

if is_true "$RUN_ATAQV_QC"; then
    printf 'status\tinterpretation\nENABLED\texperimental_ATAC_derived_not_a_CUT_pass_fail_metric\n' \
        > "$qc_root/experimental_ATAC_derived_ataqv/NOTICE.tsv"
    json_files=()
    while IFS=$'\t' read -r \
        sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
        is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
        technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
        [[ "$sample_key" == "sample_key" || "$is_control" == "TRUE" ]] && continue
        case "$genome" in hg38) organism=human ;; mm39) organism=mouse ;; *) warn "ataqv unsupported organism for $genome"; continue ;; esac
        tss="$(optional_reference_value TSS_BED "$genome")"
        [[ -n "$tss" ]] || tss="$qc_root/tss_signal_profile/reference/${genome}.tss.bed"
        [[ -s "$tss" ]] || { gtf="$(reference_value GTF "$genome")"; python3 "${SCRIPT_DIR}/prepare_tss_bed.py" "$gtf" "$tss"; }
        autosomes="$qc_root/experimental_ATAC_derived_ataqv/${genome}.autosomes.txt"
        if [[ ! -s "$autosomes" ]]; then
            samtools idxstats "$(analysis_bam_path "$sample_key")" | awk '$1!="*" {name=$1; sub(/^chr/,"",name); if(name~/^[0-9]+$/)print $1}' > "$autosomes"
        fi
        peak="${OUTPUT_DIR}/05_peaks/per_sample/${sample_key}/${primary_caller}/${sample_key}.${primary_caller}.${primary_class}.bed"
        json="$qc_root/experimental_ATAC_derived_ataqv/${sample_key}.ataqv.json.gz"
        peak_metadata="${OUTPUT_DIR}/05_peaks/per_sample/${sample_key}/peakcall_metadata.tsv"
        peak_status="$(awk -F '\t' 'NR==2 {print $5}' "$peak_metadata" 2>/dev/null || true)"
        if [[ "$peak_status" != "SUCCESS" ]]; then
            printf '{"status":"SKIPPED","reason":"primary peak status %s"}\n' "${peak_status:-MISSING}" \
                > "$qc_root/experimental_ATAC_derived_ataqv/${sample_key}.SKIPPED.json"
            warn "ataqv skipped for $sample_key because primary peak status is ${peak_status:-MISSING}"
            continue
        fi
        if ! ataqv --threads "$THREADS_ATAQV" --name "$sample_key" --metrics-file "$json" --tss-file "$tss" \
            --tss-extension "$ATAQV_TSS_EXTENSION" --autosomal-reference-file "$autosomes" --ignore-read-groups \
            --peak-file "$peak" --excluded-region-file "$blacklist" "$organism" "$(analysis_bam_path "$sample_key")" \
            >"${OUTPUT_DIR}/logs/qc/${sample_key}.ataqv.log" 2>&1; then
            printf '{"status":"FAILED","reason":"ataqv command failed"}\n' \
                > "$qc_root/experimental_ATAC_derived_ataqv/${sample_key}.FAILED.json"
            warn "ataqv failed for $sample_key"
            continue
        fi
        json_files+=("$json")
    done < "$SAMPLE_MANIFEST"
    if is_true "$GENERATE_ATAQV_VIEWER" && (( ${#json_files[@]} > 0 )); then
        mkarv -f -r calculate -p calculate "$qc_root/experimental_ATAC_derived_ataqv/viewer" "${json_files[@]}"
    fi
fi
printf 'status\tthreshold_mode\nSUCCESS\tdescriptive\n' > "$qc_root/stage_status.tsv"

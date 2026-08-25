#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(DESeq2)
    library(GenomicAlignments)
    library(GenomicRanges)
    library(Rsamtools)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6L) {
    stop("Usage: consensus_track_factors.R <sample_manifest.tsv> <cohort_id> <consensus.bed> <output_root> <table_dir> <analysis|permissive|intermediate|stringent>")
}
manifest_file <- args[[1]]; cohort_id <- args[[2]]; consensus_file <- args[[3]]
output_root <- args[[4]]; table_dir <- args[[5]]; policy <- args[[6]]
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)

metadata <- read.delim(manifest_file, stringsAsFactors=FALSE, check.names=FALSE)
metadata <- metadata[metadata$cohort_id == cohort_id & metadata$is_control == "FALSE", , drop=FALSE]
if (nrow(metadata) < 2L) stop("cohort requires at least two target biological samples")
guard <- c("genome", "assay_profile", "factor", "antibody_id", "layout", "target_class",
           "analysis_duplicate_policy", "primary_peak_caller", "primary_peak_class")
mixed <- guard[vapply(metadata[guard], function(x) length(unique(x)) != 1L, logical(1))]
if (length(mixed)) stop("mixed cohort metadata: ", paste(mixed, collapse=", "))

bed <- read.delim(consensus_file, header=FALSE, stringsAsFactors=FALSE)
if (ncol(bed) < 3L || !nrow(bed)) stop("empty/invalid consensus BED")
regions <- GRanges(bed[[1]], IRanges(as.integer(bed[[2]]) + 1L, as.integer(bed[[3]])))
region_ids <- if (ncol(bed) >= 4L) bed[[4]] else paste0(cohort_id, ".region", seq_len(nrow(bed)))
names(regions) <- region_ids

bam_path <- function(key) {
    if (policy == "analysis") return(file.path(output_root, "03_alignment", "analysis", paste0(key, ".host.analysis.bam")))
    branch <- switch(policy,
        permissive="q0_dup-retained",
        intermediate="q0_dup-removed",
        stringent="q30_dup-removed",
        stop("invalid policy: ", policy)
    )
    suffix <- gsub("_", ".", branch)
    file.path(output_root, "03_alignment", "filtered", branch, paste0(metadata$sample_key[match(key, metadata$sample_key)], ".host.", suffix, ".bam"))
}
bams <- vapply(metadata$sample_key, bam_path, character(1))
if (any(!file.exists(bams))) stop("missing policy BAMs: ", paste(bams[!file.exists(bams)], collapse=", "))
names(bams) <- metadata$sample_key

flag <- scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE,
                    isNotPassingQualityControls=FALSE, isSupplementaryAlignment=FALSE)
param <- ScanBamParam(flag=flag)
paired <- metadata$layout[[1]] == "PE"
counted <- summarizeOverlaps(
    features=regions, reads=BamFileList(bams), mode="Union",
    singleEnd=!paired, ignore.strand=TRUE, fragments=FALSE,
    inter.feature=TRUE, param=param
)
counts_matrix <- assay(counted)
colnames(counts_matrix) <- metadata$sample_key
rownames(counts_matrix) <- region_ids
if (anyNA(counts_matrix) || any(counts_matrix < 0) || any(counts_matrix != round(counts_matrix))) {
    stop("counts are not complete non-negative integers")
}
counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, , drop=FALSE]
sample_totals <- colSums(counts_matrix)
count_diagnostics <- data.frame(
    sample_key=names(sample_totals), cohort_id=cohort_id, policy=policy,
    consensus_regions_with_any_count=nrow(counts_matrix),
    consensus_count_sum=as.numeric(sample_totals),
    status=ifelse(sample_totals > 0, "NONZERO", "ZERO"),
    stringsAsFactors=FALSE
)
write.table(count_diagnostics, file.path(table_dir, "consensus_count_sums.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
if (!nrow(counts_matrix)) {
    stop("no consensus region contains counts in any cohort sample")
}
zero_samples <- names(sample_totals)[sample_totals <= 0]
if (length(zero_samples)) {
    stop("zero consensus counts for samples: ", paste(zero_samples, collapse=", "))
}
storage.mode(counts_matrix) <- "integer"
dds <- DESeqDataSetFromMatrix(counts_matrix, DataFrame(row.names=colnames(counts_matrix)), design=~1)
dds <- estimateSizeFactors(dds, type="poscounts")
sf <- sizeFactors(dds)
totals <- colSums(counts(dds))
geometric_mean <- exp(mean(log(totals)))
effective <- sf * geometric_mean
robust_scale <- 1e6 / effective
if (!isTRUE(all.equal(unname(fpm(dds, robust=TRUE)),
    unname(sweep(as.matrix(counts(dds)), 2, robust_scale, "*")),
    tolerance=1e-10, check.attributes=FALSE))) stop("robust FPM identity failed")

factor_table <- data.frame(
    sample_key=names(sf), cohort_id=cohort_id, policy=policy,
    signal_unit=ifelse(paired, "fragment", "read"), size_factor=as.numeric(sf),
    deseq2_consensus_scale=as.numeric(1/sf), consensus_count_sum=as.numeric(totals),
    cohort_geometric_mean_column_sum=geometric_mean,
    robust_effective_library_size=as.numeric(effective),
    deseq2_robust_cpm_scale=as.numeric(robust_scale), stringsAsFactors=FALSE
)
write.table(factor_table, file.path(table_dir, "normalization_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
raw <- data.frame(region_id=rownames(counts_matrix), counts_matrix, check.names=FALSE)
normalized <- data.frame(region_id=rownames(counts(dds, normalized=TRUE)), counts(dds, normalized=TRUE), check.names=FALSE)
raw_connection <- gzfile(file.path(table_dir, "raw_counts.tsv.gz"), "wt")
write.table(raw, raw_connection, sep="\t", quote=FALSE, row.names=FALSE); close(raw_connection)
norm_connection <- gzfile(file.path(table_dir, "normalized_counts.tsv.gz"), "wt")
write.table(normalized, norm_connection, sep="\t", quote=FALSE, row.names=FALSE); close(norm_connection)
writeLines(capture.output(sessionInfo()), file.path(table_dir, "session_info.txt"))

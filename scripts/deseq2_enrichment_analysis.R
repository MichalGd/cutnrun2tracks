#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 11L) stop(paste(
    "Usage: deseq2_enrichment_analysis.R <manifest> <cohort> <raw_counts.tsv.gz> <outdir>",
    "<min_replicates> <alpha> <min_abs_log2fc> <block_columns_csv> <condition_order_csv> <reference_condition> <spike_scaling.tsv|.>"
))
manifest_file <- args[[1]]; cohort_id <- args[[2]]; counts_file <- args[[3]]; outdir <- args[[4]]
min_replicates <- as.integer(args[[5]]); alpha <- as.numeric(args[[6]]); min_lfc <- as.numeric(args[[7]])
block_text <- args[[8]]; order_text <- args[[9]]; reference_condition <- args[[10]]; spike_file <- args[[11]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

metadata <- read.delim(manifest_file, stringsAsFactors=FALSE, check.names=FALSE)
metadata <- metadata[metadata$cohort_id == cohort_id & metadata$is_control == "FALSE", , drop=FALSE]
if (!nrow(metadata)) stop("cohort absent from manifest")
guard <- c("genome","assay_profile","factor","antibody_id","layout","target_class",
           "analysis_duplicate_policy","primary_peak_caller","primary_peak_class")
mixed <- guard[vapply(metadata[guard], function(x) length(unique(x)) != 1L, logical(1))]
if (length(mixed)) stop("mixed target cohort: ", paste(mixed, collapse=", "))
raw_table <- read.delim(gzfile(counts_file), stringsAsFactors=FALSE, check.names=FALSE)
region_id <- raw_table[[1]]; raw <- as.matrix(raw_table[-1]); rownames(raw) <- region_id
missing <- setdiff(metadata$sample_key, colnames(raw)); if (length(missing)) stop("missing count columns: ", paste(missing, collapse=", "))
raw <- raw[, metadata$sample_key, drop=FALSE]
if (anyNA(raw) || any(raw < 0) || any(raw != round(raw))) stop("DESeq2 input must be non-negative integer counts")
storage.mode(raw) <- "integer"

condition_counts <- table(metadata$condition)
eligible <- names(condition_counts[condition_counts >= min_replicates])
if (length(eligible) < 2L) {
    writeLines('{"status":"SKIPPED","reason":"fewer than two conditions meet replicate minimum"}', file.path(outdir,"SKIPPED.json"))
    quit(save="no", status=0)
}
observed <- unique(metadata$condition)
if (order_text != "." && nzchar(order_text)) {
    requested <- trimws(strsplit(order_text, ",", fixed=TRUE)[[1]])
    if (length(setdiff(eligible, requested))) stop("condition order omits eligible conditions")
    order <- requested[requested %in% eligible]
} else order <- observed[observed %in% eligible]
if (reference_condition != "." && nzchar(reference_condition)) {
    if (!reference_condition %in% eligible) stop("reference condition is not eligible")
    order <- c(reference_condition, setdiff(order, reference_condition))
}
metadata <- metadata[metadata$condition %in% eligible, , drop=FALSE]
metadata$condition <- factor(metadata$condition, levels=order)
raw <- raw[, metadata$sample_key, drop=FALSE]
raw <- raw[rowSums(raw) > 0, , drop=FALSE]

blocks <- character()
if (block_text != "." && nzchar(block_text)) blocks <- trimws(strsplit(block_text, ",", fixed=TRUE)[[1]])
for (block in blocks) {
    if (!block %in% names(metadata)) stop("unknown block column: ", block)
    if (any(metadata[[block]] == ".") || length(unique(metadata[[block]])) < 2L) stop("invalid block column: ", block)
    metadata[[block]] <- factor(metadata[[block]])
}
formula_text <- paste("~", paste(c(blocks, "condition"), collapse=" + "))
design_formula <- as.formula(formula_text)
rownames(metadata) <- metadata$sample_key
design_matrix <- model.matrix(design_formula, metadata)
if (qr(design_matrix)$rank < ncol(design_matrix)) stop("differential design is not full rank")

dds <- DESeqDataSetFromMatrix(raw, metadata, design_formula)
if (spike_file != ".") {
    spike <- read.delim(spike_file, stringsAsFactors=FALSE)
    spike <- spike[match(colnames(dds), spike$sample_key), , drop=FALSE]
    if (anyNA(spike$spike_observations) || any(spike$status == "FAILED")) stop("invalid spike factors for cohort")
    effective <- spike$spike_observations / spike$spikein_to_host_ratio
    sizeFactors(dds) <- effective / exp(mean(log(effective)))
    normalization_method <- "external_spikein_observations_over_declared_ratio"
} else {
    dds <- estimateSizeFactors(dds, type="poscounts")
    normalization_method <- "DESeq2_poscounts"
}
dds <- DESeq(dds, quiet=TRUE)

write.table(data.frame(sample_key=colnames(dds), condition=metadata[colnames(dds),"condition"], size_factor=sizeFactors(dds),
                       normalization_method=normalization_method), file.path(outdir,"size_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
norm <- data.frame(region_id=rownames(dds), counts(dds, normalized=TRUE), check.names=FALSE)
connection <- gzfile(file.path(outdir,"normalized_counts.tsv.gz"), "wt"); write.table(norm, connection, sep="\t", quote=FALSE, row.names=FALSE); close(connection)
file.copy(counts_file, file.path(outdir,"raw_counts.tsv.gz"), overwrite=TRUE)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
write.table(cbind(sample_key=rownames(pca), pca), file.path(outdir,"pca.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
png(file.path(outdir,"pca.png"), width=1000, height=800, res=120); print(plotPCA(vsd, intgroup="condition") + ggtitle(paste(cohort_id,"DESeq2 enrichment"))); dev.off()
png(file.path(outdir,"dispersion.png"), width=900, height=800, res=120); plotDispEsts(dds); dev.off()

comparison_rows <- list(); comparison_dir <- file.path(outdir,"comparisons"); dir.create(comparison_dir, recursive=TRUE)
index <- 0L
for (reference_index in seq_len(length(order)-1L)) for (numerator_index in (reference_index+1L):length(order)) {
    reference <- order[[reference_index]]; numerator <- order[[numerator_index]]; index <- index + 1L
    comparison_id <- paste0(numerator,"_vs_",reference); destination <- file.path(comparison_dir, comparison_id)
    dir.create(destination, recursive=TRUE)
    result <- results(dds, contrast=c("condition", numerator, reference), alpha=alpha)
    result_table <- data.frame(region_id=rownames(result), as.data.frame(result), check.names=FALSE)
    significant <- !is.na(result_table$padj) & result_table$padj <= alpha & !is.na(result_table$log2FoldChange) & abs(result_table$log2FoldChange) >= min_lfc
    result_table$direction <- ifelse(result_table$log2FoldChange > 0, paste0("higher_in_",numerator), ifelse(result_table$log2FoldChange < 0,paste0("higher_in_",reference),"unchanged"))
    all_connection <- gzfile(file.path(destination,"all.tsv.gz"),"wt"); write.table(result_table,all_connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA"); close(all_connection)
    sig_connection <- gzfile(file.path(destination,"significant.tsv.gz"),"wt"); write.table(result_table[significant,,drop=FALSE],sig_connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA"); close(sig_connection)
    comparison_rows[[index]] <- data.frame(comparison_id=comparison_id,numerator=numerator,reference=reference,tested=nrow(result_table),
        significant=sum(significant),higher_in_numerator=sum(significant & result_table$log2FoldChange>0),
        higher_in_reference=sum(significant & result_table$log2FoldChange<0),status="SUCCESS")
}
write.table(do.call(rbind,comparison_rows),file.path(outdir,"comparison_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(paste("cohort:",cohort_id),paste("design:",formula_text),paste("normalization:",normalization_method),
             paste("comparisons:",length(comparison_rows))),file.path(outdir,"summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"session_info.txt"))
saveRDS(dds,file.path(outdir,"deseq2_object.rds"))

#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DiffBind); library(rtracklayer)})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9L) stop("Usage: diffbind_analysis.R <manifest> <cohort> <consensus> <output_root> <outdir> <min_members> <alpha> <block_columns_csv> <subtract_control:true|false>")
manifest_file <- args[[1]]; cohort_id <- args[[2]]; consensus_file <- args[[3]]; output_root <- args[[4]]
outdir <- args[[5]]; min_members <- as.integer(args[[6]]); alpha <- as.numeric(args[[7]])
block_text <- args[[8]]; subtract_control <- tolower(args[[9]]) == "true"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
metadata <- read.delim(manifest_file, stringsAsFactors=FALSE, check.names=FALSE)
targets <- metadata[metadata$cohort_id == cohort_id & metadata$is_control == "FALSE", , drop=FALSE]
if (!nrow(targets)) stop("cohort absent")
condition_counts <- table(targets$condition)
eligible <- names(condition_counts[condition_counts >= min_members])
if (length(eligible) < 2L) {writeLines('{"status":"SKIPPED","reason":"insufficient replicated conditions"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
targets <- targets[targets$condition %in% eligible,,drop=FALSE]
controls <- setNames(metadata$sample_key, metadata$sample_key)
bam_path <- function(key) file.path(output_root,"03_alignment","analysis",paste0(key,".host.analysis.bam"))
sheet <- data.frame(SampleID=targets$sample_key,Tissue=targets$cell_type,Factor=targets$factor,
    Condition=targets$condition,Treatment=targets$treatment,Replicate=as.integer(targets$replicate),
    bamReads=vapply(targets$sample_key,bam_path,character(1)),
    bamControl=vapply(targets$control_key,bam_path,character(1)),ControlID=targets$control_key,
    Peaks=consensus_file,PeakCaller="bed",stringsAsFactors=FALSE)
write.csv(sheet,file.path(outdir,"diffbind_samplesheet.csv"),row.names=FALSE,quote=TRUE)
db <- dba(sampleSheet=sheet)
db <- dba.blacklist(db, blacklist=FALSE, greylist=TRUE)
db <- dba.count(db, peaks=import(consensus_file), bSubControl=subtract_control,
                bScaleControl=subtract_control, bUseSummarizeOverlaps=TRUE)
blocks <- character(); if(block_text != "." && nzchar(block_text)) blocks <- trimws(strsplit(block_text,",",fixed=TRUE)[[1]])
if(length(blocks)) {writeLines('{"status":"SKIPPED","reason":"DiffBind v0.1 does not map arbitrary block columns; use primary DESeq2Enrichment"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
design <- "~Condition"
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=min_members)
db <- dba.analyze(db, method=DBA_DESEQ2)
contrast_count <- nrow(dba.show(db,bContrasts=TRUE))
for (index in seq_len(contrast_count)) {
    result <- dba.report(db, contrast=index, method=DBA_DESEQ2, th=1)
    write.table(as.data.frame(result),file.path(outdir,paste0("contrast_",index,"_all.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    significant <- dba.report(db, contrast=index, method=DBA_DESEQ2, th=alpha)
    write.table(as.data.frame(significant),file.path(outdir,paste0("contrast_",index,"_significant.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
saveRDS(db,file.path(outdir,"diffbind_object.rds"))
writeLines(c(paste("status: SUCCESS"),paste("design:",design),paste("subtract_control:",subtract_control),paste("contrasts:",contrast_count)),file.path(outdir,"summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"session_info.txt"))

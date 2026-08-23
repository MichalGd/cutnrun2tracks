#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DESeq2);library(GenomicAlignments);library(GenomicRanges);library(Rsamtools)})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=7L)stop("Usage: control_interaction_analysis.R <manifest> <cohort> <consensus> <output_root> <outdir> <min_replicates> <alpha>")
manifest_file<-args[[1]];cohort_id<-args[[2]];consensus_file<-args[[3]];output_root<-args[[4]];outdir<-args[[5]]
minimum<-as.integer(args[[6]]);alpha<-as.numeric(args[[7]]);dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
metadata<-read.delim(manifest_file,stringsAsFactors=FALSE,check.names=FALSE)
targets<-metadata[metadata$cohort_id==cohort_id & metadata$is_control=="FALSE",,drop=FALSE]
if(length(unique(targets$condition))!=2L){writeLines('{"status":"SKIPPED","reason":"v0.1 interaction analysis requires exactly two conditions"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
if(any(table(targets$condition)<minimum) || any(targets$control_key==".")){writeLines('{"status":"SKIPPED","reason":"controls or target replicates incomplete"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
if(anyDuplicated(targets$control_key)){writeLines('{"status":"SKIPPED","reason":"interaction model requires one-to-one replicate-matched controls"}',file.path(outdir,"SKIPPED.json"));quit(save="no",status=0)}
controls<-metadata[match(targets$control_key,metadata$sample_key),,drop=FALSE]
if(anyNA(controls$sample_key) || any(controls$condition!=targets$condition)){stop("controls are not condition matched")}
bed<-read.delim(consensus_file,header=FALSE,stringsAsFactors=FALSE);regions<-GRanges(bed[[1]],IRanges(bed[[2]]+1L,bed[[3]]));names(regions)<-bed[[4]]
combined<-rbind(transform(targets,library_type="target"),transform(controls,library_type="control"))
combined$condition<-factor(combined$condition,levels=unique(targets$condition));combined$library_type<-factor(combined$library_type,levels=c("control","target"))
bams<-file.path(output_root,"03_alignment","analysis",paste0(combined$sample_key,".host.analysis.bam"));names(bams)<-paste(combined$sample_key,combined$library_type,sep=".")
if(any(!file.exists(bams)))stop("interaction BAM missing")
flag<-scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,isNotPassingQualityControls=FALSE,isSupplementaryAlignment=FALSE)
counted<-summarizeOverlaps(regions,BamFileList(bams),mode="Union",singleEnd=combined$layout[[1]]!="PE",ignore.strand=TRUE,fragments=FALSE,inter.feature=TRUE,param=ScanBamParam(flag=flag))
raw<-assay(counted);raw<-raw[rowSums(raw)>0,,drop=FALSE];storage.mode(raw)<-"integer"
rownames(combined)<-colnames(raw);design_formula<-~library_type+condition+library_type:condition
if(qr(model.matrix(design_formula,combined))$rank<ncol(model.matrix(design_formula,combined)))stop("interaction design not full rank")
dds<-DESeqDataSetFromMatrix(raw,combined,design_formula);dds<-estimateSizeFactors(dds,type="poscounts");dds<-DESeq(dds,quiet=TRUE)
coefficient<-grep("library_typetarget\\.condition",resultsNames(dds),value=TRUE)
if(length(coefficient)!=1L)stop("unable to identify interaction coefficient: ",paste(resultsNames(dds),collapse=", "))
result<-results(dds,name=coefficient,alpha=alpha);result_table<-data.frame(region_id=rownames(result),as.data.frame(result),check.names=FALSE)
connection<-gzfile(file.path(outdir,"interaction_results_all.tsv.gz"),"wt");write.table(result_table,connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(connection)
significant<-!is.na(result_table$padj)&result_table$padj<=alpha
connection<-gzfile(file.path(outdir,"interaction_results_significant.tsv.gz"),"wt");write.table(result_table[significant,,drop=FALSE],connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(connection)
writeLines(c("status: SUCCESS",paste("coefficient:",coefficient),"interpretation: condition effect in target minus condition effect in matched control","normalization: joint DESeq2 poscounts; validate before primary use"),file.path(outdir,"summary.txt"))
saveRDS(dds,file.path(outdir,"interaction_deseq2_object.rds"));writeLines(capture.output(sessionInfo()),file.path(outdir,"session_info.txt"))

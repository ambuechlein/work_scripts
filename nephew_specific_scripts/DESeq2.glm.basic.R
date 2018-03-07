#!/usr/bin/R
library('DESeq2')
library( "gplots" )
library("RColorBrewer")
library("BiocParallel")
register(MulticoreParam(64))

## args <- commandArgs(TRUE)
## control <- args[1]

cat("Loading DESeq2 data\n")
## counts <- read.table( file = args[2], header = TRUE, row.names=1 )
## pheno <- read.table( file = args[3], header=TRUE)
## rownames(pheno) <- pheno$sample
## Load pheno data matrix to creat design.
pheno <- read.table(file="phenoData2", head=TRUE, sep="\t", stringsAsFactors=TRUE, row.names=1)

# load raw counts table
counts <- read.table( file="windowCounts.tsv", header=TRUE, row.names=1 )

cat("Creating DESeq matrix and running main function\n")
dds <- DESeqDataSetFromMatrix( countData = counts, colData = pheno, design = ~ group)
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds, parallel=TRUE)
cat("1 of 14: Calculating fold-change\n")
####
p_p.aa_vs_c <- results( dds, contrast = c("group", "p_p_african_american", "p_p_caucasian"), parallel=TRUE )
p_p.aa_vs_c<-p_p.aa_vs_c[order(p_p.aa_vs_c$padj),]
write.table(as.data.frame(p_p.aa_vs_c), file="p_p_african_american_vs_p_p_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

p_p.aa_vs_h <- results( dds, contrast = c("group", "p_p_african_american", "p_p_hispanic"), parallel=TRUE )
p_p.aa_vs_h<-p_p.aa_vs_h[order(p_p.aa_vs_h$padj),]
write.table(as.data.frame(p_p.aa_vs_h), file="p_p_african_american_vs_p_p_hispanic.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

p_p.h_vs_c <- results( dds, contrast = c("group", "p_p_hispanic", "p_p_caucasian"), parallel=TRUE )
p_p.h_vs_c<-p_p.h_vs_c[order(p_p.h_vs_c$padj),]
write.table(as.data.frame(p_p.h_vs_c), file="p_p_hispanic_vs_p_p_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
####
m_p.aa_vs_c <- results( dds, contrast = c("group", "m_p_african_american", "m_p_caucasian"), parallel=TRUE )
m_p.aa_vs_c<-m_p.aa_vs_c[order(m_p.aa_vs_c$padj),]
write.table(as.data.frame(m_p.aa_vs_c), file="m_p_african_american_vs_m_p_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

m_p.aa_vs_h <- results( dds, contrast = c("group", "m_p_african_american", "m_p_hispanic"), parallel=TRUE )
m_p.aa_vs_h<-m_p.aa_vs_h[order(m_p.aa_vs_h$padj),]
write.table(as.data.frame(m_p.aa_vs_h), file="m_p_african_american_vs_m_p_hispanic.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

m_p.h_vs_c <- results( dds, contrast = c("group", "m_p_hispanic", "m_p_caucasian"), parallel=TRUE )
m_p.h_vs_c<-m_p.h_vs_c[order(m_p.h_vs_c$padj),]
write.table(as.data.frame(m_p.h_vs_c), file="m_p_hispanic_vs_m_p_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
######################
save.image(paste("Onthophagus_DESeq2.RData", sep=""))
writeLines(capture.output(sessionInfo()), "sessionInfo.deseq2.txt")

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
#p_p.aa_vs_c <- results( dds, contrast = c("group", "p_p_african_american", "p_p_caucasian"), parallel=TRUE )
#p_p.aa_vs_c<-p_p.aa_vs_c[order(p_p.aa_vs_c$padj),]
#write.table(as.data.frame(p_p.aa_vs_c), file="p_p_african_american_vs_p_p_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

m_p.h_vs_m_p.aa.c <- results(dds, contrast=list("groupm_p_hispanic", c("groupm_p_african_american","groupm_p_caucasian")), listValues=c(1, -1/2), parallel=TRUE)
m_p.h_vs_m_p.aa.c <- m_p.h_vs_m_p.aa.c[order(m_p.h_vs_m_p.aa.c$padj),]
write.table(as.data.frame(m_p.h_vs_m_p.aa.c), file="m_p_hispanic_vs_m_p_african_american_and_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

p_p.h_vs_p_p.aa.c <- results(dds, contrast=list("groupp_p_hispanic", c("groupp_p_african_american","groupp_p_caucasian")), listValues=c(1, -1/2), parallel=TRUE)
p_p.h_vs_p_p.aa.c <- p_p.h_vs_p_p.aa.c[order(p_p.h_vs_p_p.aa.c$padj),]
write.table(as.data.frame(p_p.h_vs_p_p.aa.c), file="p_p_hispanic_vs_p_p_african_american_and_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)


m_p.aa_vs_m_p.h.c <- results(dds, contrast=list("groupm_p_african_american", c("groupm_p_hispanic","groupm_p_caucasian")), listValues=c(1, -1/2), parallel=TRUE)
m_p.aa_vs_m_p.h.c <- m_p.aa_vs_m_p.h.c[order(m_p.aa_vs_m_p.h.c$padj),]
write.table(as.data.frame(m_p.aa_vs_m_p.h.c), file="m_p_african_american_vs_m_p_hispanic_and_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

p_p.aa_vs_p_p.h.c <- results(dds, contrast=list("groupp_p_african_american", c("groupp_p_hispanic","groupp_p_caucasian")), listValues=c(1, -1/2), parallel=TRUE)
p_p.aa_vs_p_p.h.c <- p_p.aa_vs_p_p.h.c[order(p_p.aa_vs_p_p.h.c$padj),]
write.table(as.data.frame(p_p.aa_vs_p_p.h.c), file="p_p_african_american_vs_p_p_hispanic_and_caucasian.DESeq2_Results.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)


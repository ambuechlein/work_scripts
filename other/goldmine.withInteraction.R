#!/usr/bin/env R
library(methylaction2)
library("BiocParallel")
library(goldmine)

args <- commandArgs(TRUE)
file <- args[1]
# file <- "all.3HH_vs_4RARA.WithInteractions.dmrs.csv"
# file <- "all.Male_vs_Female.WithInteractions.dmrs.csv"
out <- gsub(".dmrs.csv","",file)
query <- read.csv(file)

cachedir <- "/N/dc2/projects/cgbgsf/Rnor_6.0/gbcache"
genome <- "rn6"
genes <- getGenes(genome=genome,cachedir=cachedir,geneset="ensembl")
genes.r <- getGenes(genome=genome,cachedir=cachedir,geneset="refseq")

genes <- genes[str_detect(genes$isoform.id,"ENS"),]
features <- getCpgFeatures(genome=genome,cachedir=cachedir)

gm <- goldmine(query=query,genes=genes,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)
gm.r <- goldmine(query=query,genes=genes.r,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)

final <- cbind(gm$context[,c(1:5,8,11,12,14,17,18,20,23,24,26,29,30,40:47)],gm.r$context[,40:44])

if(grepl("Male_vs_Female",out)){
  colnames(final) <- c("chr","start","end","width","strand","Main.Male_vs_Female.log2FoldChange","Main.Male_vs_Female.pvalue","Main.Male_vs_Female.padj","X4RARA.Male_vs_Female.log2FoldChange","X4RARA.Male_vs_Female.pvalue","X4RARA.Male_vs_Female.padj","X3HH.Male_vs_Female.log2FoldChange","X3HH.Male_vs_Female.pvalue","X3HH.Male_vs_Female.padj","InteractionTerm.Male_vs_Female.log2FoldChange","InteractionTerm.Male_vs_Female.pvalue","InteractionTerm.Male_vs_Female.padj","ENSEMBL.call","ENSEMBL.call_genes","ENSEMBL.overlapped_genes","ENSEMBL.nearest_genes","ENSEMBL.distance_to_nearest_gene","cpgIsland_per","cpgShore_per","cpgShelf_per","Refseq.call","Refseq.call_genes","Refseq.overlapped_genes","Refseq.nearest_genes","Refseq.distance_to_nearest_gene")
}else{
  colnames(final) <- c("chr","start","end","width","strand","Main.3HH_vs_4RARA.log2FoldChange","Main.3HH_vs_4RARA.pvalue","Main.3HH_vs_4RARA.padj","Female.3HH_vs_4RARA.log2FoldChange","Female.3HH_vs_4RARA.pvalue","Female.3HH_vs_4RARA.padj","Male.3HH_vs_4RARA.log2FoldChange","Male.3HH_vs_4RARA.pvalue","Male.all.3HH_vs_4RARA.padj","InteractionTerm.3HH_vs_4RARA.log2FoldChange","InteractionTerm.3HH_vs_4RARA.pvalue","InteractionTerm.3HH_vs_4RARA.padj","ENSEMBL.call","ENSEMBL.call_genes","ENSEMBL.overlapped_genes","ENSEMBL.nearest_genes","ENSEMBL.distance_to_nearest_gene","cpgIsland_per","cpgShore_per","cpgShelf_per","Refseq.call","Refseq.call_genes","Refseq.overlapped_genes","Refseq.nearest_genes","Refseq.distance_to_nearest_gene")
}

write.table(final, file=paste(out,".annotatedResults.tsv",sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)


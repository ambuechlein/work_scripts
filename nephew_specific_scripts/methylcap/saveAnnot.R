library("BSgenome")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

anno.mart.gene.tss = MEDIPS.getAnnotation(host="www.ensembl.org", dataset=c("hsapiens_gene_ensembl"), annotation=c("TSS","GENE"), tssSz=c(-2000,500))
anno.mart.gene.tss = MEDIPS.getAnnotation(host="grch37.ensembl.org", dataset=c("hsapiens_gene_ensembl"), annotation=c("TSS","GENE"), tssSz=c(-2000,500))

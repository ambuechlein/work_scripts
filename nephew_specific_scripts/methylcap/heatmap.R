library( "gplots" )
library("RColorBrewer")
library("BSgenome")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
lnames <- load("medipsFinal.Rdata")
mr.edgeR.s = MEDIPS.selectSig(results=mr.edgeR, p.value=0.05, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
select <- order(mr.edgeR.s$edgeR.adj.p.value)
mr.edgeR.s <- mr.edgeR.s[select,]

rpkms <- mr.edgeR.s[grep(".bam.rpkm",colnames(mr.edgeR.s))]
rownames(rpkms) <- paste(mr.edgeR.s$chr,mr.edgeR.s$start,mr.edgeR.s$stop, sep="_")
colnames(rpkms) <- colN
rpkms <- rpkms[1:100,]

df <- as.data.frame(colN2)
colnames(df) <- "Condition"
rownames(df) <- colN
ann_colors <- list(Condition = c(HA="black", LA="yellow"))
pdf("tepper_heatmap.top100.pdf", height=12, width=12)
ph <- pheatmap(rpkms, color = colorRampPalette(c("navy", "white", "red"))(100), annotation_col = df, annotation_colors = ann_colors, border_color = NA, show_rownames=FALSE)
ph <- pheatmap(rpkms, color = greenred(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, show_rownames=FALSE)

heatmap.2(as.matrix(rpkms), col=greenred(256), breaks=seq(-2,2,length.out=257), dendrogram="both", scale="row", margin=c(6, 8), key=F, density.info="none", trace="none", xlab="Sample ID", cexCol=1, cexRow=.6, na.rm=TRUE,  main="Top 100 Genes", Rowv=FALSE, labRow=FALSE  )

heatmap.2(as.matrix(rpkms), scale="row", trace="none", dendrogram="both", col = greenred(256), cexCol=1, cexRow=.6, margin=c(6, 8), keysize=.8, labRow=FALSE)
heatmap.2(as.matrix(rpkms), scale="row", trace="none", dendrogram="both", col = greenred(256), cexCol=1, cexRow=.6, margin=c(6, 8), keysize=.8, labRow=FALSE, ColSideColors=colN3)
heatmap.2(as.matrix(rpkms), scale="row", trace="none", dendrogram="both", col = greenred(256), cexCol=1, cexRow=.6, margin=c(6, 8), keysize=.8, labRow=FALSE, ColSideColors=colN3, breaks=seq(-2,2,length.out=257))

dev.off()


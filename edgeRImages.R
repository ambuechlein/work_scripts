#!/usr/bin/R
args <- commandArgs(TRUE)
arg1 <- unlist(strsplit(args[1], ","))
cond1 <- arg1[1]
cond2 <- arg1[2]
group <- unlist(strsplit(args[2], ","))
outname <- paste(cond2,"_vs_",cond1, sep="")
.f = function() {
cond1 <- "MDA-MB-231_ns-tRFP_vector_FLgRNA"
cond2 <- "MDA-MB-231_MIR215_FLgRNA"
group <- unlist(strsplit("MCF7_Z_vector_FLgRNA,MCF7_Z_vector_FLgRNA,MCF7_Zeb1_FLgRNA,MCF7_Zeb1_FLgRNA,MCF7_Zeb2_FLgRNA,MCF7_Zeb2_FLgRNA,MDA-MB-231_MIR215_FLgRNA,MDA-MB-231_MIR215_FLgRNA,MDA-MB-231_ns-tRFP_vector_FLgRNA,MDA-MB-231_ns-tRFP_vector_FLgRNA,MDA231_Mir200b_FLgRNA,MDA231_Mir200b_FLgRNA,MDA231_Mir335_FLgRNA,MDA231_Mir335_FLgRNA", ","))
outname <- paste(cond2,"_vs_",cond1, sep="")
}
library(edgeR)
raw.data <- read.table( file = "countstable.txt", header = TRUE )
counts <- raw.data[ , -c(1,ncol(raw.data)+1) ]
rownames( counts ) <- raw.data[ , 1 ]
cds <- DGEList( counts , group = group )
cds <- calcNormFactors( cds )
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds , prior.df = 10 )
de.aut <- exactTest( cds, pair = c( cond1 , cond2 )) # default dispersion="auto"
de.tgw <- exactTest( cds, pair = c( cond1 , cond2 ), dispersion="tagwise")
resultsTbl.aut <- topTags( de.aut, n = nrow( de.aut$table ) )$table
de.genes.aut <- rownames( resultsTbl.aut )[ resultsTbl.aut$FDR <= 0.05 ]
summary(de <- decideTestsDGE(de.aut, p=0.05, adjust="BH"))
detags <- rownames(cds)[as.logical(de)]
resultsTbl.aut$FC <- 2^resultsTbl.aut$logFC
resultsTbl.aut <- resultsTbl.aut[c("FC", "logFC", "logCPM", "PValue", "FDR")]
wh.rows.aut <- match( rownames( resultsTbl.aut ), rownames( cds$counts ) )
wh.rows.aut2 <- match( rownames( resultsTbl.aut ), rownames( cds$pseudo.counts ) )
combResults.aut <- cbind(resultsTbl.aut, "Dispersion"=cds$tagwise.dispersion[ wh.rows.aut ], "UpDown"=decideTestsDGE( de.aut, p.value = 0.05 )[ wh.rows.aut ], cds$counts[ wh.rows.aut , ], cds$pseudo.counts[ wh.rows.aut2 , ] )
sigResultsTbl.aut <- resultsTbl.aut[ resultsTbl.aut$FDR < 0.05, ]
library(gplots)
y <- cpm(cds, normalized.lib.sizes=TRUE)
y <- y[order(match(rownames(y), rownames(resultsTbl.aut))),]
dists <- dist( t( y ) )
sigResults <- resultsTbl.aut[ resultsTbl.aut$PValue < 0.05, ]
select <- order(sigResults$PValue)
select <- seq(1, 100, 1)
pdf(paste(outname,"_images.edgeR.pdf", sep=""), pointsize=10, height=12, width=10);
heatmap.2(y[select,], col=greenred(256), dendrogram="none", symm=FALSE, scale="row", margin=c(14, 9), key=T, keysize=0.5, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=.7, na.rm=TRUE, lmat=rbind( c(0,3,4), c(2,1, 0) ), lhei=c(.5,4.5), lwid=c(.8,2,.8 ), Rowv=FALSE, Colv=FALSE )
#heatmap.2(as.matrix(dists), symm=TRUE, col=greenred(256), key=T, keysize=1, density.info="none", trace="none", labRow=cds$samples$group, margin=c(14, 14), Rowv=TRUE, Colv=TRUE )
heatmap.2(as.matrix(dists), symm=TRUE, col=greenred(256), key=T, keysize=1, density.info="none", trace="none", labRow=cds$samples$group, margin=c(14, 0), Rowv=TRUE, Colv=TRUE, cexCol=.7, cexRow=.7, lmat=rbind( c(0,3,4), c(2,1, 0) ), lhei=c(.8,4.8), lwid=c(.8,4,.8 ) )
plotSmear(de.aut, de.tags=de.genes.aut, cex = .1, ylab="Log Fold-Change", xlab="Log Counts Per Million")
abline(h = c(-2, 2), col = "blue")
library(ggplot2)
upReg <- sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == 1)
downReg <- sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == -1)
test <- data.frame("Up or Down Regulated"=factor(c("Up Regulated","Down Regulated"), levels=c("Up Regulated","Down Regulated")), Count=c(sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == 1), sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == -1)) )
ggplot(data=test, aes(x=Up.or.Down.Regulated, y=Count, fill=Up.or.Down.Regulated)) + geom_bar(stat="identity", position=position_dodge(), size=.3, width=.25) + xlab("Up or Down Regulated") + scale_fill_manual(name="Up or Down Regulated", values=c("firebrick","green4")) + theme_bw(base_size=8)
ty <- t(y)
pca.res <- prcomp(ty, center=TRUE, retx=TRUE)
plot(pca.res$x[,1:2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
legend("bottomright", row.names(ty), cex=1, ncol=1, col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
dev.off()

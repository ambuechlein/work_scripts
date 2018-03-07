#!/usr/bin/R
library(edgeR)
library(gplots)
library(ggplot2)

args <- commandArgs(TRUE)
.f = function(){
  args <- c("Mesenchymal", "/nfs/labs/nephew/SGI_Core/vs_all_tcga/Mesenchymal.c1d1_tcga_counts.tsv", "/nfs/labs/nephew/SGI_Core/vs_all_tcga/phenoDataM")
  args <- c("FLAG-vector", "/nfs/labs/nephew/nextseq_150514/fold_change/nextseq_150514.countstable.tsv", "/nfs/labs/nephew/nextseq_150514/fold_change/phenoData")
}
control <- args[1]

cat("Loading edgeR data\n")
countsTable <- read.table( file = args[2], header = TRUE, row.names=1, stringsAsFactors=TRUE )
phenoData <- read.table( file = args[3], header=TRUE)
group <- as.character(unlist(phenoData$condition))
conds <- factor(phenoData$condition)
keep <- rowSums(cpm(countsTable)>1) >= length(phenoData$condition)/2
countsTable <- countsTable[keep,]


cat("Estimating Dispersion EdgeR\n")
cds <- DGEList( countsTable, group = group )
#cds <- cds[! cds$all.zeros,]
cds <- calcNormFactors( cds )
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds , prior.df = 10 )
#cds <- estimateTagwiseDisp( cds , prior.df = 25 ) # 50/(# of samples - # of groups)
for (cond2 in unique(phenoData$condition)) {
    if(control==cond2) {
      next
    }
    cat("Creating output for ", control, "_vs_", cond2, "\n")
    outname <- paste(cond2, "_vs_", control, sep="")

    de.aut <- exactTest( cds, pair = c( control , cond2 )) # default dispersion="auto"
    de.tgw <- exactTest( cds, pair = c( control , cond2 ), dispersion="tagwise")
    resultsTbl.aut <- topTags( de.aut, n = nrow( de.aut$table ) )$table
    de.genes.aut <- rownames( resultsTbl.aut )[ resultsTbl.aut$FDR <= 0.05 ]
    summary(de <- decideTestsDGE(de.aut, p=0.05, adjust="BH"))
    detags <- rownames(cds)[as.logical(de)]
    resultsTbl.aut$FC <- 2^resultsTbl.aut$logFC
    resultsTbl.aut <- resultsTbl.aut[c("FC", "logFC", "logCPM", "PValue", "FDR")]
    write.table(resultsTbl.aut, file=paste(outname,".edgeR_Results.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
    wh.rows.aut <- match( rownames( resultsTbl.aut ), rownames( cds$counts ) )
    wh.rows.aut2 <- match( rownames( resultsTbl.aut ), rownames( cds$pseudo.counts ) )
    combResults.aut <- cbind(resultsTbl.aut, "Dispersion"=cds$tagwise.dispersion[ wh.rows.aut ], "UpDown"=decideTestsDGE( de.aut, p.value = 0.05 )[ wh.rows.aut ], cds$counts[ wh.rows.aut , ], cds$pseudo.counts[ wh.rows.aut2 , ] )
    write.table(combResults.aut, file=paste(outname,".edgeR_ResultsDetailed.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
    write.table(cds$pseudo.counts, file=paste(outname,".edgeR_NormalizedCounts.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
    sigResultsTbl.aut <- resultsTbl.aut[ resultsTbl.aut$FDR < 0.05, ]
    write.table(sigResultsTbl.aut, file=paste(outname,".edgeR_SigResults.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
    
    y <- cpm(cds, normalized.lib.sizes=TRUE)
    y <- y[order(match(rownames(y), rownames(resultsTbl.aut))),]
    dists <- dist( t( y ) )
    sigResults <- resultsTbl.aut[ resultsTbl.aut$PValue < 0.05, ]
    select <- order(sigResults$PValue)
    pdf(paste(outname,".edgeRImages.pdf", sep=""))

#    FDRvals <- c(.000001, .00001, .0001, .001, .01, .05)
#    for(i in FDRvals){
#      sigResults <- resultsTbl.aut[ resultsTbl.aut$FDR < i, ]
#      select <- order(sigResults$FDR)
#      if( length(select) >= 2){
#        heatmap.2(y[select,], col=greenred(256), dendrogram="column", symm=FALSE, scale="row", margin=c(12, 2), key=T, keysize=0.5, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=1.2, na.rm=TRUE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(1,4.5), lwid=c(1.25,2 ), labRow=FALSE, main=paste("FDR value: ", i, sep="")  )
#      }
#    }
#   
#    Pvals <- c(.000001, .00001, .0001, .001, .01, .05)
#    for(i in Pvals){
#      sigResults <- resultsTbl.aut[ resultsTbl.aut$PValue < i, ]
#      select <- order(sigResults$PValue)
#      if( length(select) >= 2){
#        heatmap.2(y[select,], col=greenred(256), dendrogram="column", symm=FALSE, scale="row", margin=c(12, 2), key=T, keysize=0.5, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=1.2, na.rm=TRUE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(1,4.5), lwid=c(1.25,2 ), labRow=FALSE, main=paste("P-value: ", i, sep="")  )
#      }
#    }

    heatmap.2(as.matrix(dists), symm=TRUE, col=greenred(256), key=T, keysize=1, density.info="none", trace="none", labRow=cds$samples$group, margin=c(14, 14), Rowv=TRUE, Colv=TRUE, main="Sample vs Sample" )

    plotSmear(de.aut, de.tags=de.genes.aut, cex = .1, ylab="Log Fold-Change", xlab="Log Counts Per Million")
    abline(h = c(-2, 2), col = "blue")

    sigResults <- resultsTbl.aut[ resultsTbl.aut$PValue < .05, ]
    select <- order(sigResults$PValue)
    if( length(select) >= 2){
      upReg <- sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == 1)
      downReg <- sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == -1)
      test <- data.frame("Up or Down Regulated"=factor(c("Up Regulated","Down Regulated"), levels=c("Up Regulated","Down Regulated")), Count=c(sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == 1), sum(decideTestsDGE(de.aut, p=0.05, adjust="BH") == -1)) )
      ggplot(data=test, aes(x=Up.or.Down.Regulated, y=Count, fill=Up.or.Down.Regulated)) + geom_bar(stat="identity", position=position_dodge(), size=.3, width=.25) + xlab("Up or Down Regulated") + scale_fill_manual(name="Up or Down Regulated", values=c("firebrick","green4")) + theme_bw(base_size=8)
    }
    Xsq <- cor(y)
#    pairs(Xsq)
    
    ty <- t(y)
    pca.res <- prcomp(ty, center=TRUE, retx=TRUE)
    screeplot(pca.res)
#    pairs(pca.res$x, col=ty)
    #biplot(pca.res, scale=0)
    plot(pca.res$x[,1:2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
    legend("topleft", row.names(ty), cex=.75, ncol=1, col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
    dev.off()
}
save.image(paste(args[2], "_edgeR.RData", sep=""))
writeLines(capture.output(sessionInfo()), "sessionInfo.edgeR.txt")

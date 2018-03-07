#!/usr/bin/R
library(DESeq)
library(gplots)

args <- commandArgs(TRUE)
# args <- c("Mesenchymal", "/nfs/labs/nephew/SGI_Core/vs_all_tcga/Mesenchymal.c1d1_tcga_counts.tsv", "/nfs/labs/nephew/SGI_Core/vs_all_tcga/phenoDataM")
# args <- c("Aa-RML12-Mock", "/nfs/projects/solexa/analysis/drosophila.hardy/foldChange/aedesCounts.tsv", "/nfs/projects/solexa/analysis/drosophila.hardy/foldChange/phenoDataA")
control <- args[1]

cat("Loading DESeq data\n")
countsTable <- read.table( file = args[2], header = TRUE, row.names=1, stringsAsFactors=TRUE )
phenoData <- read.table( file = args[3], header=TRUE)
conds <- factor(phenoData$condition)

cat("Normalizing DESeq data\n")
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
normcds <- round(as.data.frame(counts( cds, normalized=TRUE )))
write.table(normcds, file=paste(args[2],".DESeqnormalized.countstable.txt",sep=""), quote=FALSE, sep="\t", row.names=TRUE)


cat("Estimating Dispersion DESeq data\n")
# sharingMode="maximum", method="pooled" are the defaults for estimateDispersions
cds <- estimateDispersions( cds, sharingMode="maximum", method="pooled")
# For analysis without replicates
# cds <- estimateDispersions( cds, sharingMode="fit-only", method="blind")

for (cond2 in unique(phenoData$condition)) {
    if(control==cond2) {
      next
    }
    cat("Creating output for ", control, "_vs_", cond2, "\n")
    outname <- paste(cond2, "_vs_", control, sep="")

    res <- nbinomTest( cds, control, cond2 )
    res <- na.omit(res)
    
    deglist <- res[ order(res$pval), ]
    write.table(deglist, file=paste(outname,".DESeq_Results.tsv",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
    
    # 5% FDR
    resSig <- res[ res$padj < 0.05, ]
    resSig <- resSig[ order(resSig$pval), ]
    write.table(resSig, file=paste(outname,".DESeq_sig_results.tsv",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
    
    pdf(paste(outname,".DESeq.pdf",sep=""))
    
    plotDispEsts <- function( cds )
    {
        plot(
            rowMeans( counts( cds, normalized=TRUE ) ),
            fitInfo(cds)$perGeneDispEsts,
            pch = '.', log="xy", main=paste("Fit of Dispersion Estimate: ",control," vs ",cond2,sep=""), 
            xlab="Mean Expression Strength", ylab="Dispersion Values" )
        xg <- 10^seq( -.5, 5, length.out=300 )
        lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
    
    }
    plotDispEsts( cds )
    
    plotDE <- function( res )
        plot(
            res$baseMean,
            res$log2FoldChange,
            log="x", pch=20, cex=.3,
            main= paste("Log2 Fold Change Versus Base Means (5% FDR): ",control," vs ",cond2,sep=""), sub="MA-Plot",
            xlab="Base Means", ylab="Log2 Fold Change",
            col = ifelse( res$padj < .05, "red", "black" ) )
    plotDE( res )
    
    hist(res$pval, breaks=100, col="skyblue", border="slateblue", main=paste("Histogram of p-values: ",control," vs ",cond2,sep=""), xlab="p-value")
    
    cdsBlind <- estimateDispersions( cds, method="blind" )
    vsd <- getVarianceStabilizedData( cdsBlind )
    
    mod_lfc <- (rowMeans( vsd[, conditions(cds)==control, drop=FALSE] ) - rowMeans( vsd[, conditions(cds)==cond2, drop=FALSE] ))
    lfc <- res$log2FoldChange
    finite <- is.finite(lfc)
    table(as.character(lfc[!finite]), useNA="always")
    largeNumber <- 10
    lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)
    
    logdecade <- 1 + round( log10( 1+rowMeans(counts(cdsBlind, normalized=TRUE)) ) )
#    colors <- colorRampPalette( c( "gray", "blue" ) )(6)[logdecade]
#    plot( lfc, mod_lfc, pch=20, cex=.4, asp=1,
#          main=paste("Scatterplot of direct (lfc) versus moderated log-ratios (mod_lfc): ",control," vs ",cond2,sep=""), 
#          col = ifelse( finite, colors, "purple" ) )
#    abline( a=0, b=1, col="#40C04040" )
    
    dists <- dist( t( vsd ) )
    heatmap.2(as.matrix(dists), symm=TRUE, col=greenred(256), key=T, keysize=1, density.info="none", trace="none", labRow = paste(pData(cdsBlind)$condition, pData(cdsBlind)$type), margin=c(12, 12), Rowv=TRUE, Colv=TRUE, main="Sample vs Sample" )
    
    select <- order(res$pval)[1:40]
    heatmap.2( vsd[select,], col = greenred(256), scale = "none", main="Top 40 Most Significant Genes")
    
    ##########################################
     FDRvals <- c(.000001, .00001, .0001, .001, .01, .05)
     for(i in FDRvals){
       sigResults <- res[ res$padj < i, ]
       select <- order(sigResults$padj)
       if( length(select) >= 2){
         heatmap.2(vsd[select,], col=greenred(256), dendrogram="column", symm=FALSE, scale="row", margin=c(12, 2), key=T, keysize=0.5, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=1.2, na.rm=TRUE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(1,4.5), lwid=c(1.25,2 ), labRow=FALSE, main=paste("FDR value: ", i, sep="")  )
       }
     }
    
     Pvals <- c(.000001, .00001, .0001, .001, .01, .05)
     for(i in Pvals){
       sigResults <- res[ res$pval < i, ]
       select <- order(sigResults$pval)
       if( length(select) >= 2){
         heatmap.2(vsd[select,], col=greenred(256), dendrogram="column", symm=FALSE, scale="row", margin=c(12, 2), key=T, keysize=0.5, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=1.2, na.rm=TRUE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(1,4.5), lwid=c(1.25,2 ), labRow=FALSE, main=paste("P-value: ", i, sep="")  )
       }
     }
    ##########################################
    
    Xsq <- cor(vsd)
    pairs(Xsq)
    
    tvsd <- t(vsd)
    pca.res <- prcomp(tvsd, center=TRUE, retx=TRUE)
    plot(pca.res$x[,1:2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
    legend("topright", row.names(tvsd), cex=.6, ncol=2, col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), pch=c(1,2,3,4,5,6,7,8,9,10,11,12,15,17,19))
    dev.off()
}
save.image(paste(args[2], "_DESEeq.RData", sep=""))
writeLines(capture.output(sessionInfo()), "sessionInfo.deseq.txt")

#!/usr/bin/R
library('DESeq2')
library( "gplots" )
library("RColorBrewer")
library("BiocParallel")
library(pheatmap)
register(MulticoreParam(32))

args <- commandArgs(TRUE)
# args <- c("Norm", "Ivan_Hypoxia.countstable.tsv", "phenoData")
control <- args[1]

cat("Loading DESeq2 data\n")
counts <- read.table( file = args[2], header = TRUE, row.names=1 )
phenoData <- read.table( file = args[3], header=TRUE, sep="\t", stringsAsFactors=TRUE,)
rownames(phenoData) <- phenoData$sample

cat("Creating DESeq matrix and running main function\n")
dds <- DESeqDataSetFromMatrix( countData = counts, colData = phenoData, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 5, ] 
#dds <- dds[ rowSums(counts(dds)[,grep("C1D1",colnames(dds))])/15 > 5 | rowSums(counts(dds)[,grep("C2D8",colnames(dds))])/15 > 5, ]
dds$condition <- relevel( dds$condition, control )

dds <- DESeq(dds, parallel=TRUE)
for (cond2 in unique(phenoData$condition)) {
    if(control==cond2) {
      next
    }
    cat("Creating output for ", control, "_vs_", cond2, "\n")
    outname <- paste(cond2, "_vs_", control, sep="")

    pdf(paste(outname, ".deseq2.pdf", sep=""))
cat("1 of 14: Calculating fold-change\n")
    res <- results( dds, contrast = c("condition", cond2, control), parallel=TRUE )
    select <- order(res$padj)
    select05 <- na.pass(res$padj < 0.05)
    select05 <- select05 %in% TRUE
    filter <- na.pass(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
    filter <- filter %in% TRUE
    sigResults <- na.omit(res[ na.omit(res$padj) < .05, ])
    resFull <- results( dds, cooksCutoff=FALSE, independentFiltering=FALSE, contrast = c("condition", cond2, control), parallel=TRUE )
    res.row <- match( rownames( res ), rownames( resFull ) )
    resAll <- cbind(as.data.frame(res), "UnfilteredFDR"=resFull$padj[ res.row ])
    colnames(resAll) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "FilteredFDR", "UnfilteredFDR")
#    res<-res[order(res$padj),]
#    resFull <- resFull[order(resFull$padj),]
    write.table(as.data.frame(res)[select,], file=paste(outname,".DESeq2_Results.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
    write.table(as.data.frame(resFull)[order(resFull$padj),], file=paste(outname,".DESeq2_Results.full.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
    write.table(resAll, file=paste(outname,".DESeq2_Results.all.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)

cat("1 of 14: Creating MAplot\n")
    plotMA( res, ylim = c(-1, 1) )
    hist( res$pvalue, breaks=20, col="grey" )
    
    # create bins using the quantile function
cat("3 of 14: create bins using the quantile function\n")
    qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
    # "cut" the genes into the bins
    bins <- cut( res$baseMean, qs )
    # rename the levels of the bins using the middle point
    levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
    # calculate the ratio of $p$ values less than .01 for each bin
    ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
    # plot these ratios
cat("4 of 14: plotting ratio barplot\n")
    barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
cat("5 of 14: plotting basemean quantiles\n")
    plot(metadata(res)$filterNumRej,type="b", ylab="number of rejections", xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)

    #library( "biomaRt" )
    #ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
    #genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
    #                  filters = "ensembl_gene_id",
    #                  values = rownames(res),
    #                  mart = ensembl )
    #idx <- match( rownames(res), genemap$ensembl_gene_id )
    #res$entrez <- genemap$entrezgene[ idx ]
    #res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

  cat("6 of 14: plotting Dispersion Estimates\n")
      plotDispEsts( dds, ylim = c(1e-6, 1e1) )
  cat("7 of 14: calculating rlog\n")
      rld <- rlog( dds )
      vsd <- varianceStabilizingTransformation(dds)
      rlogMat <- assay(rld)
      vstMat <- assay(vsd)
#      library("vsn")
#      par(mfrow=c(1,3))
#      notAllZero <- (rowSums(counts(dds))>0)
#      meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
#      meanSdPlot(rlogMat[notAllZero,])
#      meanSdPlot(vstMat[notAllZero,])
      par( mfrow = c( 1, 2 ) )
  cat("8 of 14: plotting log2\n")
      plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
  cat("9 of 14: plotting rlog\n")
      plot( rlogMat[, 1:2], col="#00000020", pch=20, cex=0.3 )
  cat("10 of 14: Calculating Distance Matrix\n")
      sampleDists <- dist( t( rlogMat ) )
      sampleDistMatrix <- as.matrix( sampleDists )
      rownames(sampleDistMatrix) <- rld$sample
      colnames(sampleDistMatrix) <- rld$sample
  cat("11 of 14: Plotting sample vs sample heatmap\n")
      heatmap.2(sampleDistMatrix, symm=TRUE, col=greenred(256), key=T, keysize=1, density.info="none", trace="none", margin=c(14, 14), Rowv=TRUE, Colv=TRUE, main="Sample vs Sample" )
  cat("12 of 14: Plotting PCA\n")
      print( plotPCA( rld, intgroup = c( "condition", "sample")))
      data <- plotPCA(rld, intgroup=c("condition", "sample"), returnData=TRUE)
      percentVar <- round(100 * attr(data, "percentVar"))
      print(ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")))
      print(ggplot(data, aes(PC1, PC2, color=sample, shape=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")))

  cat("13 of 14: finding top variable genes\n")
      topVarGenes <- head( order( genefilter::rowVars( rlogMat ), decreasing=TRUE ), 35 )
  cat("14 of 14: plotting heatmap for top variable genes\n")
      heatmap.2( rlogMat[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
      
      select <- order(res$padj)
      resTMP <- res[select,]
      vstPlot <- vstMat[select,]
      rlogPlot <- rlogMat[select,]
      select05 <- na.pass(resTMP$padj < 0.05)
      select05 <- select05 %in% TRUE
      if(nrow(vstPlot[select05,]) > 500){
        select05[501:length(select05)] <- FALSE
      } 
      heatmap.2( vstPlot[select05,], col = greenred(256), dendrogram="none", symm=FALSE, scale="none", Rowv=FALSE, margin=c(6,6), key=F, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=.7, na.rm=TRUE, main="FDR < 0.05: Variance Stabalized Values", Colv=FALSE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(.2,1), lwid=c(1,6), cexRow=.7 )
      heatmap.2( rlogPlot[select05,], col = greenred(256), dendrogram="none", symm=FALSE, scale="none", Rowv=FALSE, margin=c(6,6), key=F, density.info="none", trace="none", xlab="Sample ID", ylab= "Gene ID", cexCol=.7, na.rm=TRUE, main="FDR < 0.05: R-Log Values", Colv=FALSE, lmat=rbind( c(4,3), c(2,1) ), lhei=c(.2,1), lwid=c(1,6), cexRow=.7 )

      mat <- rlogMat[select[1:50],]
      mat <- mat - rowMeans(mat)
      colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
      colnames(mat) <- phenoData$sample
      df <- as.data.frame(colData(rld)[,c("condition")])
#      df <- as.data.frame(colnames(mat))
      colnames(df) <- "CellType"
      rownames(df) <- colnames(mat)
      df2 <- gsub("\\d+_","", df$CellType)
      df2 <- gsub("-","_",df2)
#      df$CellType <- df2
      #ann_colors <- list(CellType = c(C1D1="blue", C2D8="red"))
      mycols <- brewer.pal(8, "Dark2")[1:length(unique(phenoData$condition))]
      ann_colors <- list(CellType=c(mycols))
      names(ann_colors$CellType) <- unique(phenoData$condition)
      
      title <- "Top 50 Significant Heatmap"
      ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row", cluster_rows=FALSE)
      ###

      mat <- rlogMat[filter,]
      mat <- mat - rowMeans(mat)
      colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
      colnames(mat) <- phenoData$sample
      df <- as.data.frame(colData(rld)[,c("condition")])
#      df <- as.data.frame(colnames(mat))
      colnames(df) <- "CellType"
      rownames(df) <- colnames(mat)
      df2 <- gsub("\\d+","", df$CellType)
      df2 <- gsub("-","_",df2)
#      df$CellType <- df2
      #ann_colors <- list(CellType = c(LC="blue", LD="red", TB="orange", TC="yellow", TD="green", TT="black"))
      title <- "FDR < 0.05 and FC > 2 Heatmap"
      ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row", cluster_rows=FALSE)
      ###

      mat <- rlogMat[ topVarGenes, ]
      mat <- mat - rowMeans(mat)
      colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
      colnames(mat) <- phenoData$sample
      df <- as.data.frame(colData(rld)[,c("condition")])
#      df <- as.data.frame(colnames(mat))
      colnames(df) <- "CellType"
      rownames(df) <- colnames(mat)
      df2 <- gsub("\\d+_","", df$CellType)
      df2 <- gsub("-","_",df2)
#      df$CellType <- df2
      #ann_colors <- list(CellType = c(C1D1="blue", C2D8="red"))
      title <- "Top Variable Genes Heatmap"
      ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row", cluster_rows=FALSE)
    dev.off()
    cat(outname, " output complete\n")
}
sess_out <- gsub(".tsv", "",args[2])
sess_out <- gsub(".countstable", "", sess_out)
save.image(paste(sess_out, "_DESeq2.RData", sep=""))
writeLines(capture.output(sessionInfo()), paste(sess_out,".sessionInfo.deseq2.txt", sep=""))

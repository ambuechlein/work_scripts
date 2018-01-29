#!/usr/bin/R
library('DESeq2')
library( "gplots" )
library("ggplot2")
library("RColorBrewer")
library("BiocParallel")
library(pheatmap)
library("IHW")
register(MulticoreParam(32))

args <- commandArgs(TRUE)
control <- args[1]

cat("Loading DESeq2 data\n")
counts <- read.table( file = args[2], header = TRUE, row.names=1 )
phenoData <- read.table( file = args[3], header=TRUE, sep="\t", stringsAsFactors=TRUE,)
rownames(phenoData) <- phenoData$sample

cat("Creating DESeq matrix and running main function\n")
dds <- DESeqDataSetFromMatrix( countData = counts, colData = phenoData, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 5, ] 
#dds <- dds[ rowSums(counts(dds)[,grep("C1D1",colnames(dds))])/15 > 5 | rowSums(counts(dds)[,grep("C2D8",colnames(dds))])/15 > 5, ]
dds$condition <- relevel( dds$condition, ref=control )

dds <- DESeq(dds, parallel=TRUE)
cat("Calculating rlog and vst\n")
rld <- rlog( dds )
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

sess_out <- gsub(".tsv", "",args[2])
sess_out <- gsub(".countstable", "", sess_out)
write.table(counts(dds, normalized=TRUE),file=paste(sess_out,".normalizedCounts.tsv",sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
for (cond2 in unique(phenoData$condition)) {
  if(control==cond2) {
    next
  }
  cat("Creating output for ", control, "_vs_", cond2, "\n")
  outname <- paste(cond2, "_vs_", control, sep="")

  pdf(paste(outname, ".deseq2.pdf", sep=""))
  cat("1 of 14: Calculating fold-change\n")
  res <- results( dds, contrast = c("condition", cond2, control), parallel=TRUE )
#### IHW ####
  resDF <- as.data.frame(res)
  ihw_res <- ihw(pvalue ~ baseMean,  data=resDF, alpha = 0.1)
  ihw_res_df <- as.data.frame(ihw_res)
  fullResDF <- cbind(resDF,ihw_res_df)
  selectIHW <- order(fullResDF$adj_pvalue)
  filterIHW <- na.pass(fullResDF$adj_pvalue < 0.05 & abs(fullResDF$log2FoldChange) > 1)
  filterIHW <- filterIHW %in% TRUE
  write.table(fullResDF[selectIHW,], file=paste(outname,".DESeq2_Results.IHW.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
  plotDF <- fullResDF
  plotDF$color <- "black"
  plotDF$color[plotDF$adj_pvalue < 0.05] = "blue"
  plotDF$color[plotDF$padj < 0.05] = "green"
  plotDF$color[plotDF$adj_pvalue < 0.05 & plotDF$padj < 0.05] = "red"
  pdf(paste(outname,".IHW_plots.pdf",sep=""))
  plot(plotDF$adj_pvalue,plotDF$padj, ylim=c(0,1), col=plotDF$color, xlim=c(0,1), main="Standard vs IHW P-value Correction", xlab="IHW Adjusted P-value", ylab="DESeq2 Default Adjusted P-value")
  legend("topleft", legend=c("Not Significant","IHW Significant","Default Significant","Both Significant"), col=c("black","blue","green","red"), pch=20)
  
  plot(plotDF$adj_pvalue,plotDF$pvalue, ylim=c(0,1), col=plotDF$color, xlim=c(0,1), main="IHW vs P-value", xlab="IHW Adjusted P-value", ylab="P-value")
  legend("topleft", legend=c("Not Significant","IHW Significant","Default Significant","Both Significant"), col=c("black","blue","green","red"), pch=20)
  
  plot(plotDF$padj,plotDF$pvalue, ylim=c(0,1), col=plotDF$color, xlim=c(0,1), main="Standard vs P-value", xlab="DESeq2 Adjusted P-value", ylab="P-value")
  legend("topleft", legend=c("Not Significant","IHW Significant","Default Significant","Both Significant"), col=c("black","blue","green","red"), pch=20)
  dev.off()
#############

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
  write.table(as.data.frame(res)[select,], file=paste(outname,".DESeq2_Results.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
  write.table(as.data.frame(resFull)[order(resFull$padj),], file=paste(outname,".DESeq2_Results.full.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
  write.table(resAll, file=paste(outname,".DESeq2_Results.all.tsv", sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)

  # create bins using the quantile function
  cat("2 of 14: create bins using the quantile function\n")
  qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
  # "cut" the genes into the bins
  bins <- cut( res$baseMean, qs )
  # rename the levels of the bins using the middle point
  levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
  # calculate the ratio of $p$ values less than .01 for each bin
  ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
  # plot these ratios
  cat("3 of 14: plotting ratio barplot\n")
  barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
  cat("4 of 14: plotting basemean quantiles\n")
  plot(metadata(res)$filterNumRej,type="b", ylab="number of rejections", xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)

  #library( "biomaRt" )
  #ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
  #genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
  #          filters = "ensembl_gene_id",
  #          values = rownames(res),
  #          mart = ensembl )
  #idx <- match( rownames(res), genemap$ensembl_gene_id )
  #res$entrez <- genemap$entrezgene[ idx ]
  #res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

  cat("5 of 14: plotting Dispersion Estimates\n")
  plotDispEsts( dds, ylim = c(1e-6, 1e1) )
  par( mfrow = c( 1, 2 ) )
  cat("6 of 14: plotting log2\n")
  plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
  cat("7 of 14: plotting rlog\n")
  plot( rlogMat[, 1:2], col="#00000020", pch=20, cex=0.3 )
  cat("8 of 14: Creating MAplot\n")
  par( mfrow = c( 1, 1 ) )
  hist( res$pvalue, breaks=20, col="grey" )
#  plotMA( res, ylim = c(-1, 1) )
  plotMA( res, ylim = c(-1, 1), alpha=0.05)
  cat("9 of 14: Creating volcano plot\n")
  volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
    with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
    with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
    if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
    }
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
  }
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
#  volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
#  volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2),labelsig=FALSE)
  volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, labelsig=FALSE)

resdata$threshold[abs(resdata$log2FoldChange) > 1] <- "|LogFC| > 1"
resdata$threshold[resdata$padj < 0.05] <- "FDR < 0.05"
resdata$threshold[abs(resdata$log2FoldChange) > 1 & resdata$padj < 0.05] <- "Both"
resdata$threshold <- factor(resdata$threshold,levels=c("Not Significant","|LogFC| > 1","FDR < 0.05","Both"))

print(ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot All"))

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
  shapes <- as.integer(unique(phenoData$condition))
  if(length(shapes) > 6){
    print(ggplot(data, aes(PC1, PC2, color=sample, shape=condition)) + scale_shape_manual(values=shapes) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")))
  } else {
    print(ggplot(data, aes(PC1, PC2, color=sample, shape=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")))
  }

  cat("13 of 14: finding top variable genes\n")
  topVarGenes <- head( order( genefilter::rowVars( rlogMat ), decreasing=TRUE ), 35 )
  cat("14 of 14: plotting heatmap for top variable genes\n")
  select <- order(res$padj)
  resTMP <- res[select,]
  vstPlot <- vstMat[select,]
  rlogPlot <- rlogMat[select,]
  select05 <- na.pass(resTMP$padj < 0.05)
  select05 <- select05 %in% TRUE
  mat <- rlogMat[select[1:50],]
  mat <- mat - rowMeans(mat)
  colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
  colnames(mat) <- phenoData$sample
  df <- as.data.frame(colData(rld)[,c("condition")])
  colnames(df) <- "CellType"
  rownames(df) <- colnames(mat)
  df2 <- gsub("\\d+_","", df$CellType)
  df2 <- gsub("-","_",df2)
  mycols <- brewer.pal(8, "Dark2")[1:length(unique(phenoData$condition))]
  ann_colors <- list(CellType=c(mycols))
  names(ann_colors$CellType) <- unique(phenoData$condition)

  title <- "Top 50 Significant Heatmap"
  if(nrow(mat) > 2){
    ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
    ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
  }
###
  mat <- rlogMat[filter,]
  if( !is.null(nrow(mat)) && nrow(mat) > 2){
    mat <- mat - rowMeans(mat)
    colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
    colnames(mat) <- phenoData$sample
    df <- as.data.frame(colData(rld)[,c("condition")])
    colnames(df) <- "CellType"
    rownames(df) <- colnames(mat)
    df2 <- gsub("\\d+","", df$CellType)
    df2 <- gsub("-","_",df2)
    title <- "FDR < 0.05 and FC > 2 Heatmap"
    if( !is.null(nrow(mat)) && nrow(mat) > 2){
      ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
      ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
    }
  }
###
  mat <- rlogMat[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  colnames(mat) <- paste(phenoData$condition, " ", phenoData$sample, sep="")
  colnames(mat) <- phenoData$sample
  df <- as.data.frame(colData(rld)[,c("condition")])
  colnames(df) <- "CellType"
  rownames(df) <- colnames(mat)
  df2 <- gsub("\\d+_","", df$CellType)
  df2 <- gsub("-","_",df2)
  title <- "Top Variable Genes Heatmap"
  if( !is.null(nrow(mat)) && nrow(mat) > 2){
    ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
    ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=FALSE, cluster_cols=FALSE, scale="row")
  }
  dev.off()
  cat(outname, " output complete\n")
}
save.image(paste(sess_out, "_DESeq2.RData", sep=""))
writeLines(capture.output(sessionInfo()), paste(sess_out,".sessionInfo.deseq2.txt", sep=""))

library("BiocParallel")
register(MulticoreParam(32))
library(DESeq2)
library(sva)
library(devtools)
library(limma)
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
library(VennDiagram)
library(plyr)
library(dplyr)
library(tidyr)
library(methylaction2)
library(goldmine)

lnames <- load("all.3HH_vs_4RARA.WithInteractions.RData")
rm(reads,windows,bins,test.one,wingr,bins.gr,filter.pass,zero)
gc()
resForPlotting <- as.data.frame(ma$dmr)

resForPlotting$threshold <- "Not Significant"
resForPlotting$threshold[abs(resForPlotting$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange) > 1] <- "|LogFC| > 1"
resForPlotting$threshold[resForPlotting$Main.all.3HH_vs_4RARA.WithInteractions.padj < 0.05] <- "FDR < 0.05"
resForPlotting$threshold[abs(resForPlotting$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange) > 1 & resForPlotting$Main.all.3HH_vs_4RARA.WithInteractions.padj < 0.05] <- "Both"
resForPlotting$threshold <- factor(resForPlotting$threshold,levels=c("Not Significant","|LogFC| > 1","FDR < 0.05","Both"))
resForPlotting$bmRank <- rank(resForPlotting$Main.all.3HH_vs_4RARA.WithInteractions.baseMean)
var <- length(resForPlotting$bmRank)/3
resForPlotting$BaseMeanRank <- "Bottom Third"
resForPlotting$BaseMeanRank[resForPlotting$bmRank > var+var] <- "Top Third"
resForPlotting$BaseMeanRank[resForPlotting$bmRank < var+var & resForPlotting$bmRank > var] <- "Middle Third"
resForPlotting$BaseMeanRank <- factor(resForPlotting$BaseMeanRank,levels=c("Top Third","Middle Third","Bottom Third"))

rownames(resForPlotting) <- paste(resForPlotting$seqnames,":",resForPlotting$start,"..",resForPlotting$end,sep="")
colnames(resForPlotting)[1] <- "chr"
resForPlotting$chr <- paste("chr",resForPlotting$chr,sep="")

rld <- rlog( cds )
rlogMat <- assay(rld)
rownames(rlogMat) <- rownames(resForPlotting)

cachedir <- "/N/dc2/projects/cgbgsf/Rnor_6.0/gbcache"
genome <- "rn6"
genes <- getGenes(genome=genome,cachedir=cachedir,geneset="ensembl")
genes.r <- getGenes(genome=genome,cachedir=cachedir,geneset="refseq")

genes <- genes[str_detect(genes$isoform.id,"ENS"),]
features <- getCpgFeatures(genome=genome,cachedir=cachedir)

gm <- goldmine(query=resForPlotting,genes=genes,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)
gm.r <- goldmine(query=resForPlotting,genes=genes.r,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)
#final <- cbind(gm$context,gm.r$context[,40:44])
gm$context$pattern <- as.character(gm$context$pattern)
gm.r$context$pattern <- as.character(gm.r$context$pattern)

gm$context$pattern[gm$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange > 0] <- "Hypermethylated"
gm$context$pattern[gm$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange < 0] <- "Hypomethylated"

gm$context$pattern[gm$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange > 0 & gm$context$Main.all.3HH_vs_4RARA.WithInteractions.padj < .05] <- "Hypermethylated and Significant"
gm$context$pattern[gm$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange < 0 & gm$context$Main.all.3HH_vs_4RARA.WithInteractions.padj < .05] <- "Hypomethylated and Significant"

gm.r$context$pattern[gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange > 0] <- "Hypermethylated"
gm.r$context$pattern[gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange < 0] <- "Hypomethylated"

gm.r$context$pattern[gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange > 0 & gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.padj < .05] <- "Hypermethylated and Significant"
gm.r$context$pattern[gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange < 0 & gm.r$context$Main.all.3HH_vs_4RARA.WithInteractions.padj < .05] <- "Hypomethylated and Significant"

sigResults <- na.omit(subset(resForPlotting,Main.all.3HH_vs_4RARA.WithInteractions.padj < .05))

pdf(paste(out,".QCPlots.pdf",sep=""))

print( plotPCA( rld, intgroup = c( "group", "batch")))
data <- plotPCA(rld, intgroup=c("group", "batch","sex"), returnData=TRUE)
data$SA <- paste(data$group.1,":",data$batch,sep="")

percentVar <- round(100 * attr(data, "percentVar"))
shapes <- as.integer(unique(colData(rld)$group))
print(ggplot(data, aes(PC1, PC2, color=group.1,shape=sex)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")))

sampleDists <- dist( t( rlogMat ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$batch,".",rld$sample,sep="")
colnames(sampleDistMatrix) <- paste(rld$batch,".",rld$sample,sep="")
title <- "Sample vs Sample Distance Heatmap"
ph <- pheatmap(sampleDistMatrix, color = colorRampPalette(c("navy", "white", "red"))(256), border_color = NA, main=title, fontsize_row=4, fontsize_col=4,)

sampleDists <- dist( t( rlogMat[gm$context$call=="promoter",] ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$batch,".",rld$sample,sep="")
colnames(sampleDistMatrix) <- paste(rld$batch,".",rld$sample,sep="")
title <- "Sample vs Sample Distance Heatmap Promoter Only"
ph <- pheatmap(sampleDistMatrix, color = colorRampPalette(c("navy", "white", "red"))(256), border_color = NA, main=title, fontsize_row=4, fontsize_col=4,)

ggplot(data=resForPlotting, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot All")

ggplot(data=resForPlotting, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), )) +
  geom_point(alpha=0.4, size=1.75,aes(colour=log10(Main.all.3HH_vs_4RARA.WithInteractions.baseMean))) + scale_colour_gradientn(colours=rainbow(4),name="log10(BaseMean)")+
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot All") 

ggplot(data=resForPlotting, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot All")

ggplot(data=subset(resForPlotting,  BaseMeanRank == "Top Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot Top Third Basemean Rank")

ggplot(data=subset(resForPlotting,  BaseMeanRank == "Middle Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot Middle Third Basemean Rank")

ggplot(data=subset(resForPlotting,  BaseMeanRank == "Bottom Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot Bottom Third Basemean Rank")

ggplot(data=sigResults, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR < 0.05")

ggplot(data=sigResults, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR < 0.05")

ggplot(data=subset(sigResults,  BaseMeanRank == "Top Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR < 0.05 Top Third Basemean Rank")

ggplot(data=subset(sigResults,  BaseMeanRank == "Middle Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR < 0.05 Middle Third Basemean Rank")

ggplot(data=subset(sigResults,  BaseMeanRank == "Bottom Third"), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=BaseMeanRank)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR < 0.05 Bottom Third Basemean Rank")

###
gencon <- gm$context[,list(count=length(chr)),by=c("pattern","call")]
gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
gencon <- gencon[,list(call=call,count=count,total=sum(count),
                        fraction=count/sum(count)),by="pattern"]

genconHyper <- gencon[grep("Hypermethylated",gencon$pattern),]
genconHypo <- gencon[grep("Hypomethylated",gencon$pattern),]

genconHyperSig <- gencon[gencon$pattern=="Hypermethylated and Significant",]
genconHypoSig <- gencon[gencon$pattern=="Hypomethylated and Significant",]

ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of DMRs ENSEMBL")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs ENSEMBL")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs ENSEMBL") + coord_polar("y")


ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs ENSEMBL")
ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs ENSEMBL") + coord_polar("y")

ggplot(genconHyperSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypermethylated and Significant DMRs ENSEMBL")
ggplot(genconHyperSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypermethylated and Significant DMRs ENSEMBL") + coord_polar("y")

ggplot(genconHypoSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypomethylated and Significant DMRs ENSEMBL")
ggplot(genconHypoSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypomethylated and Significant DMRs ENSEMBL") + coord_polar("y")

featcon <- gm$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                         CPGshore=sum(cpgShore_per>0)/length(chr),
                         CPGshelf=sum(cpgShelf_per>0)/length(chr)),
                         by=c("pattern")]
featcon <- melt(featcon,id.vars=c("pattern"))
setnames(featcon,c("variable","value"),c("call","percent"))

ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                        position="dodge") + ggnice() + labs(title="Feature Context of DMRs ENSEMBL")

gencon <- gm.r$context[,list(count=length(chr)),by=c("pattern","call")]
gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
gencon <- gencon[,list(call=call,count=count,total=sum(count),
                        fraction=count/sum(count)),by="pattern"]

genconHyper <- gencon[grep("Hypermethylated",gencon$pattern),]
genconHypo <- gencon[grep("Hypomethylated",gencon$pattern),]

genconHyperSig <- gencon[gencon$pattern=="Hypermethylated and Significant",]
genconHypoSig <- gencon[gencon$pattern=="Hypomethylated and Significant",]

ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of DMRs Refseq")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs Refseq")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs Refseq") + coord_polar("y")


ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs Refseq")
ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs Refseq") + coord_polar("y")

ggplot(genconHyperSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypermethylated and Significant DMRs Refseq")
ggplot(genconHyperSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypermethylated and Significant DMRs Refseq") + coord_polar("y")

ggplot(genconHypoSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypomethylated and Significant DMRs Refseq")
ggplot(genconHypoSig,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypomethylated and Significant DMRs Refseq") + coord_polar("y")

featcon <- gm.r$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                         CPGshore=sum(cpgShore_per>0)/length(chr),
                         CPGshelf=sum(cpgShelf_per>0)/length(chr)),
                         by=c("pattern")]
featcon <- melt(featcon,id.vars=c("pattern"))
setnames(featcon,c("variable","value"),c("call","percent"))

ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                        position="dodge") + ggnice() + labs(title="Feature Context of DMRs Refseq")

p<-ggplot(gm$context, aes(x=Main.all.3HH_vs_4RARA.WithInteractions.baseMean, fill=call, color=call)) + geom_histogram(position="dodge", alpha=0.5) + xlim(c(0,1500)) + ylim(c(0,4000)) + labs(title="Histogram of Basemean by Genomic Type")
# Add mean lines
mu <- ddply(gm$context, "call", summarise, grp.mean=mean(Main.all.3HH_vs_4RARA.WithInteractions.baseMean))
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=call),linetype="dashed")

ggplot(subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.baseMean, fill=call, color=call)) + geom_histogram(position="dodge", alpha=0.5) + labs(title="Histogram of Basemean of Promoters Only")

ggplot(subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.baseMean, fill=call, color=call)) + geom_histogram(position="dodge", alpha=0.5) + xlim(c(0,100)) + labs(title="Histogram of Basemean of Promoters Only")

ggplot(data=subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR Promoter Only")

ggplot(data=subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.log2FoldChange, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj),)) +
  geom_point(alpha=0.4, size=1.75,aes(colour=log10(Main.all.3HH_vs_4RARA.WithInteractions.baseMean))) + scale_colour_gradientn(colours=rainbow(4),name="log10(BaseMean)") +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + geom_vline(xintercept = 1, linetype="dashed") + geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2 fold change") + ylab("-log10 fdr") + labs(title="Volcano Plot FDR Promoter Only")

ggplot(data=subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.baseMean, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + #geom_vline(xintercept = mean(Main.all.3HH_vs_4RARA.WithInteractions.baseMean), linetype="dashed") +
  xlab("BaseMean") + ylab("-log10 fdr") + labs(title="BaseMean vs -log10FDR")

ggplot(data=subset(gm$context, call=="promoter" ), aes(x=Main.all.3HH_vs_4RARA.WithInteractions.baseMean, y=-log10(Main.all.3HH_vs_4RARA.WithInteractions.padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) + xlim(c(0,75)) +
  geom_hline(yintercept = -log10(.05), linetype="dashed") + #geom_vline(xintercept = mean(Main.all.3HH_vs_4RARA.WithInteractions.baseMean), linetype="dashed") +
  xlab("BaseMean") + ylab("-log10 fdr") + labs(title="BaseMean vs -log10FDR")

###

phenoData$SA <- paste(phenoData$sex,phenoData$group,sep='.')
phenoData$SA <- factor(phenoData$SA)
allgroups <- levels(phenoData[,"SA"])
tissuecols <- brewer.pal(8, "Dark2")[1:length(unique(phenoData$SA))]
col.nut2 <- tissuecols

# Get group levels and define outliers.
# what database will be used for expression levels
nt <- normTransform(cds)
nt.assay <- assay(nt)
levelsdb <- nt.assay
# database for outliers
cooks <- assays(cds)[["cooks"]]
# generate databases with expression levels per group, and the selection database
groupmeans <- NULL
groupmaxs <- NULL
groupmins <- NULL
groupmedians <- NULL
groupranges <- NULL
maxcooks <- NULL
# define all the groups for iteration
allgroups <- levels(phenoData[,"SA"])
for(g in allgroups){
  # select all samples from this group
  selsamp <- row.names(phenoData[phenoData$SA==g,])
  # make a subdb from the levels
  subdb <- levelsdb[,selsamp]
  subcooks <- cooks[,selsamp]
  # extract variation in a row
  av <- rowMeans(subdb)
  med <- rowMedians(subdb)
  mi <- rowMin(subdb)
  mma <- rowMax(subdb)
  ran <- mma-mi
  # extract max cooks
  mc <- rowMax(subcooks)
  # put all values in respective dataframes
  groupmeans <- cbind(groupmeans,av)
  groupmaxs <- cbind(groupmaxs,mma)
  groupmins <- cbind(groupmins,mi)
  groupmedians <- cbind(groupmedians,med)
  groupranges <- cbind(groupranges,ran)
  maxcooks <- cbind(maxcooks,mc)
}
# finalize the dataframes
finalize <- function(dm){
  colnames(dm) <- allgroups
  row.names(dm) <- row.names(levelsdb)
  df <- as.data.frame(dm)
  return(df)
}

groupmeans <- finalize(groupmeans)
groupmaxs <- finalize(groupmaxs)
groupmins <- finalize(groupmins)
groupmedians <- finalize(groupmedians)
groupranges <- finalize(groupranges)
maxcooks <- finalize(maxcooks)

# Plots
# specify subsets in the groups of samples
## groupsHH <- grep("HH",allgroups,value=T)
## groupsTH <- grep("TH",allgroups,value=T)
## groupsGEN <- grep("GEN",allgroups,value=T)
## groupsBR <- grep("BR",allgroups,value=T)
groupsYes <- grep("3HH",allgroups,value=T)
groupsNo <- grep("4RARA",allgroups,value=T)
## groupsFHC <- grep("female.3HH.Coord_Blood_Sample",allgroups,value=T)
## groupsFHP <- grep("female.3HH.Peripheral_Blood",allgroups,value=T)
## groupsFLC <- grep("female.4RARA.Coord_Blood_Sample",allgroups,value=T)
## groupsFLP <- grep("female.4RARA.Peripheral_Blood",allgroups,value=T)
## groupsMHC <- grep("male.3HH.Coord_Blood_Sample",allgroups,value=T)
## groupsMHP <- grep("male.3HH.Peripheral_Blood",allgroups,value=T)
## groupsMLC <- grep("male.4RARA.Coord_Blood_Sample",allgroups,value=T)
## groupsMLP <- grep("male.4RARA.Peripheral_Blood",allgroups,value=T)
## mlist <- list(groupsFHC,groupsFHP,groupsFLC,groupsFLP,groupsMHC,groupsMHP,groupsMLC,groupsMLP)

## minimums
#plot(density(groupmins[,groupsYes[1]],na.rm=TRUE),col="white",main="minimum density")
#for (name in names(myList)) {
#
#}
plot(density(groupmins[,groupsYes[1]],na.rm=TRUE),col="white",main="minimum density")
for(g in groupsYes){
  lines(density(groupmins[,g],na.rm=TRUE),col=tissuecols[2])
}
for(g in groupsNo){
  lines(density(groupmins[,g],na.rm=TRUE),col=tissuecols[1])
}
legend("topright", legend=c("3HH","4RARA"), col=c(tissuecols[2],tissuecols[1]),pch=20)
## maximums
plot(density(groupmaxs[,groupsYes[1]],na.rm=TRUE),col="white",main="maximum density")
for(g in groupsYes){
  lines(density(groupmaxs[,g],na.rm=TRUE),col=tissuecols[2])
}
for(g in groupsNo){
  lines(density(groupmaxs[,g],na.rm=TRUE),col=tissuecols[1])
}
legend("topright", legend=c("3HH","4RARA"), col=c(tissuecols[2],tissuecols[1]),pch=20)

## means
plot(density(groupmeans[,groupsYes[1]],na.rm=TRUE),col="white",main="mean density")
for(g in groupsYes){
  lines(density(groupmeans[,g],na.rm=TRUE),col=tissuecols[2])
}
for(g in groupsNo){
  lines(density(groupmeans[,g],na.rm=TRUE),col=tissuecols[1])
}
legend("topright", legend=c("3HH","4RARA"), col=c(tissuecols[2],tissuecols[1]),pch=20)

## medians
plot(density(groupmedians[,groupsYes[1]],na.rm=TRUE),col="white",main="median density")
for(g in groupsYes){
  lines(density(groupmedians[,g],na.rm=TRUE),col=tissuecols[2])
}
for(g in groupsNo){
  lines(density(groupmedians[,g],na.rm=TRUE),col=tissuecols[1])
}
legend("topright", legend=c("3HH","4RARA"), col=c(tissuecols[2],tissuecols[1]),pch=20)

## ranges
plot(density(groupranges[,groupsYes[1]],na.rm=TRUE),col="white",main="range density")
for(g in groupsYes){
  lines(density(groupranges[,g],na.rm=TRUE),col=alpha(tissuecols[2],0.4))
}
for(g in groupsNo){
  lines(density(groupranges[,g],na.rm=TRUE),col=alpha(tissuecols[1],0.4))
}
legend("topright", legend=c("3HH","4RARA"), col=c(tissuecols[2],tissuecols[1]),pch=20)

dev.off()

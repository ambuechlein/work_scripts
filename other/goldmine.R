#!/nfs/bio/sw/bin/R
library(methylaction2)
library("BiocParallel")
library(goldmine)

args <- commandArgs(TRUE)
file <- args[1]
out <- gsub(".dmrs.csv","",file)
query <- read.csv(file)

cachedir <- "/nfs/labs/nephew/human_databases/gbcache"
genome <- "rn6"
genes <- getGenes("gencode",genome=genome,cachedir=cachedir,geneset="ensembl")
genes <- genes[str_detect(genes$isoform.id,"ENS"),]
features <- getCpgFeatures(genome=genome,cachedir=cachedir)

gm <- goldmine(query=query,genes=genes,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)

gm$context$pattern[gm$context[[grep("log2FoldChange",colnames(gm$context))]] > 0] <- "Hypermethylated"
gm$context$pattern[gm$context[[grep("log2FoldChange",colnames(gm$context))]] < 0] <- "Hypomethylated"
#gm$context$pattern[gm$context$pattern==1] <- "Hypermethylated"
#gm$context$pattern[gm$context$pattern==10] <- "Hypomethylated"

write.table(gm$context, file=paste(out,".annotatedResults.tsv",sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
write.table(gm$genes, file=paste(out,".genes.tsv",sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)
gmWrite(gm, path=paste(out,"_gm_csv",sep=""))

gencon <- gm$context[,list(count=length(chr)),by=c("pattern","call")]
gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
gencon <- gencon[,list(call=call,count=count,total=sum(count),
                        fraction=count/sum(count)),by="pattern"]

#gm$context$pattern[gm$context$pattern==1] <- "Hypermethylated"
#gm$context$pattern[gm$context$pattern==10] <- "Hypomethylated"

#gencon$pattern[gencon$pattern==1] <- "Hypermethylated"
#gencon$pattern[gencon$pattern==10] <- "Hypomethylated"

genconHyper <- gencon[gencon$pattern=="Hypermethylated",]
genconHypo <- gencon[gencon$pattern=="Hypomethylated",]
pdf(paste(out,"golminePlots.pdf",sep=""))
ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of DMRs")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs")
ggplot(genconHyper,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypermethylated DMRs") + coord_polar("y")

ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity", position="dodge") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs")
ggplot(genconHypo,aes(x="",y=fraction,fill=call)) + geom_bar(stat="identity") + ggnice() + labs(title="Genomic Context of Hypomethylated DMRs") + coord_polar("y")

featcon <- gm$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                         CPGshore=sum(cpgShore_per>0)/length(chr),
                         CPGshelf=sum(cpgShelf_per>0)/length(chr)),
                         by=c("pattern")]
featcon <- melt(featcon,id.vars=c("pattern"))
setnames(featcon,c("variable","value"),c("call","percent"))

ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                        position="dodge") + ggnice() + labs(title="Feature Context of DMRs")
dev.off()

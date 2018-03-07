library(goldmine)
query <- read.csv("dmrs.csv")
head(query)
cachedir <- "/nfs/labs/nephew/human_databases/gbcache"
genome <- "hg38"
genes <- getGenes("gencode",genome=genome,cachedir=cachedir,gencodetable="wgEncodeGencodeBasicV24")
genes <- genes[str_detect(genes$isoform.id,"ENS"),]
features <- getCpgFeatures(genome=genome,cachedir=cachedir)

gm <- goldmine(query=query,genes=genes,genome=genome,cachedir=cachedir, promoter = c(2000, 500), features=features)
write.table(gm$context, file="methyAction.Goldmine.annotatedResults.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
write.table(gm$genes, file="methyAction.Goldmine.genes.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
gmWrite(gm, path="gm_csv")

gencon <- gm$context[,list(count=length(chr)),by=c("pattern","call")]
gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
gencon <- gencon[,list(call=call,count=count,total=sum(count),
                        fraction=count/sum(count)),by="pattern"]
pdf("golminePlots.pdf")
ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity",
                        position="dodge") + ggnice() + labs(title="Genomic Context of DMRs")

featcon <- gm$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                         CPGshore=sum(cpgShore_per>0)/length(chr),
                         CPGshelf=sum(cpgShelf_per>0)/length(chr)),
                         by=c("pattern")]
featcon <- melt(featcon,id.vars=c("pattern"))
setnames(featcon,c("variable","value"),c("call","percent"))

ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                        position="dodge") + ggnice() + labs(title="Feature Context of DMRs")
dev.off()

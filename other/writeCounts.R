library("BiocParallel")
register(MulticoreParam(32))
library(DESeq2)

lnames <- load("all.3HH_vs_4RARA.WithInteractions.RData")
rm(reads,windows,bins,test.one,wingr,bins.gr,filter.pass,zero)
gc()
resForPlotting <- as.data.frame(ma$dmr)
rownames(resForPlotting) <- paste(resForPlotting$seqnames,":",resForPlotting$start,"..",resForPlotting$end,sep="")

nc <- counts(cds,normalized=FALSE)
rownames(nc) <- rownames(resForPlotting)
write.table(nc,file="RatBlood.RawCount.tsv",sep="\t",quote=FALSE)

nc <- counts(cds,normalized=TRUE)
rownames(nc) <- rownames(resForPlotting)
write.table(nc,file="RatBlood.NormalizedCount.tsv",sep="\t",quote=FALSE)


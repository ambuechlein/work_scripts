library("BSgenome")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

chrset <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
uniq <- 1e-3

lnames <- load("loadedBams.Rdata")

# get Annotation
#load("anno.mart.gene.Rdata")
load("anno.mart.gene.tss.Rdata")

#anno.mart.gene = MEDIPS.getAnnotation(dataset=c("hsapiens_gene_ensembl"), annotation=c("GENE"), )
#anno.mart.gene.tss = MEDIPS.getAnnotation(host="www.ensembl.org", dataset=c("hsapiens_gene_ensembl"), annotation=c("TSS","GENE"), tssSz=c(-2000,500))
#anno.mart.gene = MEDIPS.getAnnotation(dataset=c("hsapiens_gene_ensembl"), annotation=c("GENE"), )

# create coupling set 
# any of the other MEDIPS SETs would be fine, because all of them consist of the same set of chromosomes and have been generated with the same window size
cat("Creating Coupling Vector\n")
CS = MEDIPS.couplingVector(pattern = "CG", refObj = LAset[[1]])

### tmp for Landen to simplify
#LAset <- SsetStrict
#HAset <- LsetStrict
outname <- "tepper"
###

# differential methylation
cat("Performing ttest\n")
#mr.ttest = MEDIPS.meth(MSet1=HAset, MSet2=LAset, CSet=CS, p.adj="bonferroni", diff.method="ttest", MeDIP=T, CNV=F, type="rpkm", minRowSum=10)
cat("Performing edgeR\n")
mr.edgeR = MEDIPS.meth(MSet1=HAset, MSet2=LAset, CSet=CS, p.adj="bonferroni", diff.method="edgeR", MeDIP=T, CNV=F, type="rpkm", minRowSum=10)

# Fails because no CSet is provided
# mr.edgeR2 = MEDIPS.meth(MSet1=PSet, MSet2=RSet, p.adj="bonferroni", diff.method="edgeR", prob.method="poisson", MeDIP=T, CNV=F, type="rpkm", minRowSum=1)
# mr.edgeR3 = MEDIPS.meth(MSet1=PSet, MSet2=RSet, p.adj="bonferroni", diff.method="edgeR", prob.method="poisson", MeDIP=F, CNV=F, type="rpkm", minRowSum=10)
cat("Annotating Results\n")
#mr.ttest = MEDIPS.setAnnotation(regions=mr.ttest, annotation=anno.mart.gene.tss)
mr.edgeR = MEDIPS.setAnnotation(regions=mr.edgeR, annotation=anno.mart.gene.tss)
cat("Saving Results Rdata\n")
save.image(paste(outname, ".ranDiffMeth.Rdata", sep=""))

cat("writing results tables\n")
cat("Filtering on fdr\n")
res <- mr.edgeR
keep <- complete.cases(res$edgeR.adj.p.value)
res <- res[keep,]
res <-res [order(res$edgeR.adj.p.value),]
cat("Removing extra columns\n")
cols <- grep("^chr|start|stop|counts|rpkm|edgeR|^1_id|^2_id|^3_id", colnames(res))
res <- res[,cols]
write.table(res, file=paste(outname,".edgeR_Results.final.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
#write.table(mr.edgeR, file=paste(outname,".edgeR_Results.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)

mr.edgeR.s = MEDIPS.selectSig(results=mr.edgeR, p.value=0.05, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
if(length(mr.edgeR.s$chr) > 0){
  res.s = MEDIPS.selectSig(results=res, p.value=0.05, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
  write.table(res.s, file=paste(outname,".edgeR_sigResults.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
}

cat("select regions of interest CpG\n")
columns = names(mr.edgeR)[grep("counts|rpkm|edgeR|^1_id|^2_id|^3_id", colnames(res))]
#columns = names(mr.edgeR)[grep("counts", names(mr.edgeR))]
# select regions of interest CpG
cpg.coords <- read.table(file="/nfs/bio/db/Homo_sapien/gencode_v19/cpg_islands.forMEDIPS.bed", header=FALSE, stringsAsFactors = FALSE)
colnames(cpg.coords) <- c("chr","start","stop","ID")
cpg.rois   = MEDIPS.selectROIs(results=mr.edgeR, rois=cpg.coords, columns=columns, summarize=NULL)
cpg.rois.s = MEDIPS.selectROIs(results=mr.edgeR, rois=cpg.coords, columns=columns, summarize="avg")
#cpg.rois = MEDIPS.setAnnotation(regions=cpg.rois, annotation=anno.mart.gene.tss)
write.table(cpg.rois, file=paste(outname,".edgeR.Results.cpg.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
write.table(cpg.rois.s, file=paste(outname,".edgeR.Results.cpg.avg.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
save.image(paste(outname, ".medipsFinal.Rdata", sep=""))
MEDIPS.plotCalibrationPlot(CSet=CS, main="Calibration Plot", MSet=LAset[[1]], plot_chr="chr22", rpkm=TRUE, xrange=TRUE)
MEDIPS.plotCalibrationPlot(CSet=CS, main="Calibration Plot", MSet=HAset[[1]], plot_chr="chr22", rpkm=TRUE, xrange=TRUE)
writeLines(capture.output(sessionInfo()), paste(outname, ".sessionInfo.medips.txt",sep=""))

for(set in LAset){
  out <- sample_name(set)
  out <- gsub(".bam", ".wig", out)
  MEDIPS.exportWIG(Set=set, file=out, format="rpkm", descr="")
}
for(set in HAset){
  out <- sample_name(set)
  out <- gsub(".bam", ".wig", out)
  MEDIPS.exportWIG(Set=set, file=out, format="rpkm", descr="")
} 

if(length(mr.edgeR.s$chr) > 0){
  cat("creating various outputs\n")
  # Merge Frames/windows
  # Higher in MSet1 (Control)
  mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] > 0), ]
  mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.gain, distance=1)

  # select regions of interest merged frames (gain in Control)
  rois.g   = MEDIPS.selectROIs(results=mr.edgeR, rois=mr.edgeR.s.gain.m, columns=columns, summarize=NULL)
  rois.g.s = MEDIPS.selectROIs(results=mr.edgeR, rois=mr.edgeR.s.gain.m, columns=columns, summarize="avg")
  write.table(rois.g, file=paste(outname,".edgeR.sigResults.mergedWindows.gainInControl.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
  write.table(rois.g.s, file=paste(outname,".edgeR.sigResults.mergedWindows.gainInControl.avg.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
  
  # Merge Frames/windows
  #Lower in MSet1 (Control)
  mr.edgeR.s.loss = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] < 0), ]
  mr.edgeR.s.loss.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.loss, distance=1)
  
  # select regions of interest merged frames (loss in Control)
  rois.l   = MEDIPS.selectROIs(results=mr.edgeR, rois=mr.edgeR.s.loss.m, columns=columns, summarize=NULL)
  rois.l.s = MEDIPS.selectROIs(results=mr.edgeR, rois=mr.edgeR.s.loss.m, columns=columns, summarize="avg")
  write.table(rois.l, file=paste(outname,".edgeR.sigResults.mergedWindows.lossInControl.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
  write.table(rois.l.s, file=paste(outname,".edgeR.sigResults.mergedWindows.lossInControl.avg.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
}
save.image(paste(outname, ".medipsFinal.Rdata", sep=""))
writeLines(capture.output(sessionInfo()), paste(outname, ".sessionInfo.medips.txt",sep=""))

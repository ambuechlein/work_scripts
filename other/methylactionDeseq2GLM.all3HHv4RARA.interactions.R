#!/usr/bin/env R
library(methylaction2)
library("BiocParallel")

lnames <- load("preprocess.glm.RData")
##########################################################################################
out <- "all.3HH_vs_4RARA.WithInteractions"
#out <- "test.sex_interaction"
poifdr=0.1
stageone.p=0.05
anodev.p=0.05
post.p=0.05
freq=2/3
minsize=150
joindist=200
adjust.var=NULL
nperms=0
perm.boot=F
ncore <- 40 
register(MulticoreParam(ncore))
cov <- NULL
stagetwo.method <- "co"
winsize <- width(counts)[1]

ngroups <- length(unique(samp$group))
phenoData <- data.frame(sample=samp$sample,group=samp$group,sex=samp$gender,tissue=samp$tissue,batch=samp$batch,stringsAsFactors = TRUE)
rownames(phenoData) <- phenoData$sample
write.table(phenoData, file="p1.tsv",sep="\t", quote=FALSE,na="NA", col.names=NA)
phenoData <- read.table( file = "p1.tsv", header=TRUE, sep="\t", stringsAsFactors=TRUE,)
phenoData <- phenoData[-1]
rownames(phenoData) <- phenoData$sample
system("rm p1.tsv")

message("Starting analysis for a ",ngroups," group comparision")
samp <- data.table(samp)
if(sum(samp[,length(sample),by=group]$V1 > 1)!=ngroups){stop("Each group must contain more than one replicate")}
groups <- unique(samp$group)
if(!is.null(adjust.var)){if(!(adjust.var %in% colnames(samp))){stop("No column with name equal to adjust.var found in samp")}}

args <- list(samp=samp, stagetwo.method=stagetwo.method, winsize=winsize, poifdr=poifdr, stageone.p=stageone.p, freq=freq, joindist=joindist, anodev.p=anodev.p, adjust.var=adjust.var, post.p=post.p, minsize=minsize, nperms=nperms, perm.boot=perm.boot, ncore=ncore, start=Sys.time())

# Do Initial Filtering
fdr.filter <- methylaction2:::filter(counts, samp, poifdr)

# Get Signal Bins Only
message("Filtering to Signal Windows")
filter.pass <- rowSums(t(t(as.matrix(values(counts)))>fdr.filter$cuts))>0

# Compute Size Factors based on Filtered Regions - will always use these for all subsequent tests
message("Computing size factors")
#        sizefactors <- estimateSizeFactorsForMatrix(as.matrix(values(counts[filter.pass])))
#phenoData$SA <- paste(phenoData$sex,phenoData$group,sep='.')
dds <- DESeqDataSetFromMatrix( countData = as.matrix(values(counts[filter.pass])), colData = phenoData, design = ~ batch + sex + group)
dds$group <- relevel( dds$group, ref="4RARA")
dds$sex <- relevel( dds$sex, ref="female")

cat("Fold change\n")
dds <- DESeq(dds, parallel=TRUE)

comps <- methylaction2:::getGroupComps(unique(samp$group))

res <- results( dds, contrast = c("group","3HH","4RARA"), parallel=TRUE )
res <- list(res)
names(res) <- c(out)

sizefactors <- dds$sizeFactor
normcounts <- counts(dds,normalized=TRUE)
rm(dds)
gc()
zero <- rowSums(as.matrix(values(counts)))==0
wingr <- GRanges(seqnames(counts),IRanges(start(counts),end(counts)))
windows <- list()
windows$zero <- wingr[zero]
windows$filtered <- wingr[!filter.pass & !zero]
windows$signal <- counts[filter.pass]
windows$signal.norm <- wingr[filter.pass]
values(windows$signal.norm) <- normcounts

patt <- methylaction2:::callPatternsN(res=res,cutoff=stageone.p,samp=samp)
colgroups <- lapply(unique(samp$group),function(x) samp[group==x,]$sample)
names(colgroups) <- unique(samp$group)
counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(as.matrix(values(windows$signal.norm[,i])))))
colnames(counts.means) <- paste0(levels(samp$group),".mean")

# Combine so we get these in the output
patt <- data.table(counts.means,patt)
rm(counts,counts.means)
gc()
# join adjacent equivalent patterns
message("Reduction by pattern and disjoing regions")
bins <- windows$signal
bins.gr <- bins
values(bins.gr) <- NULL
bins.gr$pattTestOne <- patt$patt

# split out by pattern
allone <- paste(as.character(rep(1,length(unique(samp$group)))),collapse="")
sigpatt <- bins.gr[!(bins.gr$pattTestOne %in% c("ambig","000or111",allone))]
pattrows <- sigpatt$pattTestOne
values(sigpatt) <- NULL
bins.bypatt <- split(sigpatt,pattrows)

# reduce each pattern with itself within the gap distance
bins.red <- lapply(bins.bypatt,reduce,min.gapwidth=joindist)
for(i in 1:length(bins.red)){bins.red[[i]]$pattTestOne <- names(bins.red)[i]}

# only keep reduced ranges if minsize or larger
bins.red <- lapply(bins.red,function(x) x[width(x)>=minsize])

# deal with conflicting extensions using disjoin()
bins.red.gr <- do.call(c,unname(bins.red))
regions.gr <- disjoin(bins.red.gr)
regions <- regions.gr[width(regions.gr)>=minsize]

# output: genomic ranges to pass to stage 2
patterns <- bins
values(patterns) <- patt

test.one <- list(patterns=patterns, regions=regions)

###
regions=test.one$regions
stagetwo.method="co"
message("Begin stage two testing")
recounts <- as.matrix(values(getCounts(samp=samp,reads=reads,ranges=regions,ncore=ncore)))

testDESeqANODEV<- function(recounts,groups,covar, ncore)
{
cds <- DESeqDataSetFromMatrix( countData = recounts, colData = phenoData, design = ~ batch + sex + group)
sizeFactors(cds) <- sizefactors
message("Estimating dispersions")
cds <- suppressWarnings(estimateDispersions(cds,fitType="local"))
message("Fitting full GLM")
message("Fitting reduced GLM")
message("Performing test")
cds <- nbinomLRT(cds, reduced = ~ 1)
res <- results(cds,parallel=TRUE)
finalres <- as.numeric(res$pvalue)
return(finalres)
}
anodev <- testDESeqANODEV(recounts=recounts,groups=samp$group,covar=covar,ncore=ncore)
anodev.padj <- p.adjust(anodev,method="fdr")

regions$anodev.p <- anodev
regions$anodev.padj <- anodev.padj
cds <- DESeqDataSetFromMatrix( countData = recounts, colData = phenoData, design = ~ batch + sex + group)
cds$group <- relevel( cds$group, ref="4RARA")
cds$sex <- relevel( cds$sex, ref="female")
cds <- DESeq(cds, parallel=TRUE)
normcounts <- counts(cds,normalized=TRUE)
rm(cds)
gc()
anodev.keep <- (anodev.padj<anodev.p) & !(is.na(anodev.padj))
if(sum(anodev.keep)==0){message("Found no DMRs in stage two ANODEV")}
if(sum(anodev.keep)==0){out<-"NoDmrs"; return(out);}

test.two <- list()
test.two$ns <- regions[!anodev.keep]
test.two$ns.counts <- regions[!anodev.keep]
values(test.two$ns.counts) <- normcounts[!anodev.keep,]

regions.sig <- regions[anodev.keep]
recounts.sig <- recounts[anodev.keep,,drop=F]
cds <- DESeqDataSetFromMatrix( countData = recounts.sig, colData = phenoData, design = ~ batch + sex + group)
cds$group <- relevel( cds$group, ref="4RARA")
cds$sex <- relevel( cds$sex, ref="female")
cds <- DESeq(cds, parallel=TRUE)

testres <- results( cds, contrast = c("group","3HH","4RARA"), parallel=TRUE )
#testres <- list(testres)
#names(testres) <- c(out)

xds <- DESeqDataSetFromMatrix( countData = recounts.sig, colData = phenoData, design = ~ batch + sex + group + sex:group)
xds$group <- relevel( xds$group, ref="4RARA")
xds$sex <- relevel( xds$sex, ref="female")
xds <- DESeq(xds, parallel=TRUE)
# the group effect for sex female (the main effect)
f.3HHv4RARA <- results(xds, name="group_3HH_vs_4RARA",parallel=TRUE)
# the group effect for sex male
# this is, by definition, the main effect *plus* the interaction term
# (the extra group effect in sex male compared to sex female).
m.3HHv4RARA <- results(xds, list( c("group_3HH_vs_4RARA","sexmale.group3HH") ),parallel=TRUE)
# the interaction term, answering: is the group effect *different* across sexs?
i.3HHv4RARA <- results(xds, name="sexmale.group3HH",parallel=TRUE)

testres <- list(testres,f.3HHv4RARA,m.3HHv4RARA,i.3HHv4RARA)
names(testres) <- c(paste("Main.",out,sep=""), paste("Female.",out,sep=""),paste("Male.",out,sep=""),paste("InteractionTerm.",out,sep=""))

if(!all(sapply(testres,nrow)==nrow(recounts.sig))){stop("Ran out of memory during testing, try reducing ncore")}

patt <- methylaction2:::callPatternsN(res=testres, cutoff=post.p, samp=samp)

test.two$sig <- regions[anodev.keep]
values(test.two$sig) <- patt
test.two$sig.normcounts <- regions[anodev.keep]
values(test.two$sig.normcounts) <- normcounts[anodev.keep,,drop=F]

dmr <- data.frame(chr=seqnames(regions.sig),start=start(regions.sig),end=end(regions.sig),width=width(regions.sig),anodev.padj=regions.sig$anodev.padj,pattern=patt$patt)

dmrcalled <- cbind(dmr[,-5],testres)
dmrcalled$dmrid <- 1:nrow(dmrcalled)

colgroups <- lapply(unique(samp$group),function(x) samp[group==x,]$sample)
names(colgroups) <- unique(samp$group)
counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(recounts.sig[,i,drop=F])))
colnames(counts.means) <- paste0(unique(samp$group),".mean")

dmr <- cbind(dmr,counts.means,recounts.sig)

test.two$dmr <- makeGRanges(dmr)
test.two$dmrcalled <- makeGRanges(dmrcalled)
## 

ma <- list(dmr=test.two$dmrcalled, args=args, data=list(windows=windows, fdr.filter=fdr.filter, sizefactors=sizefactors, test.one=test.one, test.two=test.two))
# Remove some things from the output to reduce the size
ma$data$test.two$dmrcalled <- NULL
values(ma$data$windows$filtered) <- NULL

# Output results
ma$args$end <- Sys.time()
ma$args$hours <- as.numeric(difftime(ma$args$end,ma$args$start,units="hours"))

##########################################################################################
save.image(paste(out,".RData",sep=""), compress=TRUE)

write.csv(makeDT(ma$dmr), row.names=FALSE, file=paste(out,".dmrs.csv",sep=""))
write.table(makeDT(ma$dmr), file=paste(out,".methyAction.DMR.tsv",sep=""), quote=FALSE, sep="\t", na="NA", col.names=NA)

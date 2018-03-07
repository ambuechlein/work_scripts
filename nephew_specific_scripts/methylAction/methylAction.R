#!/nfs/bio/sw/bin/R
library(methylaction)
samp <- readSampleInfo("samples.csv")
chrs <- paste0("chr",c(1:22,"X","Y"))
fragsize <- 75
winsize <- 50
ncore <- 32
reads <- getReads(samp=samp, chrs=chrs, fragsize=fragsize, ncore=ncore)
counts <- getCounts(samp=samp, reads=reads, chrs=chrs, winsize=winsize, ncore=ncore)
save.image("preprocess.RData", compress=TRUE)
ma <- methylaction(samp=samp, counts=counts, reads=reads, ncore=ncore)
save(ma,file="ma.rd", compress=T)
ma2 <- methylaction(samp=samp, counts=counts, reads=reads, perm.boot=T, nperms=3, ncore=ncore)
save(ma2,file="ma2.rd", compress=T)
save.image("final.RData", compress=TRUE)

write.csv(makeDT(ma$dmr), row.names=FALSE, file="dmrs.csv")
write.table(makeDT(ma$dmr), file="methyAction.MA.DMR.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
write.table(makeDT(ma2$dmr), file="methyAction.MA2.DMR.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
write.table(maSummary(ma), file="methyAction.summary.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)
maSummary(ma)
maSummary(ma2)
table(ma$dmr$pattern,ma$dmr$frequent)
table(ma2$dmr$pattern,ma2$dmr$frequent)
pdf("methylAction.images.pdf")
maHeatmap(ma)
maHeatmap(ma2)
maKaryogram(ma=ma, reads=reads)
maKaryogram(ma=ma2, reads=reads)

dev.off()

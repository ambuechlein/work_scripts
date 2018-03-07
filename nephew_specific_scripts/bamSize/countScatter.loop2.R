#!/usr/bin/R

wd <- getwd()

samC.70 <- read.table( file = paste(wd,"/countsAll.tsv", sep=""), header = TRUE, row.names=1, stringsAsFactors=TRUE )
nrC.70  <- read.table( file = paste(wd,"/countsDeduped.tsv", sep=""), header = TRUE, row.names=1, stringsAsFactors=TRUE )

pomNR <- read.table( file = "/nfs/projects/solexa/analysis/human.pomerening/tophat/merged_samples_results/HAD0G2_AACGTGAT.deduped.htseq.tsv", header = TRUE, row.names=1, stringsAsFactors=TRUE )
pom <- read.table( file = "/nfs/projects/solexa/analysis/human.pomerening/tophat/merged_samples_results/HAD0G2_AACGTGAT.sorted.htseq.tsv", header = TRUE, row.names=1, stringsAsFactors=TRUE )

## plot(pom$AllReads, pomNR$Nonredundant, main="RNA-Seq", xlab="All Reads", ylab="Nonredundant Reads", pch=19)
## abline(lm(pomNR$Nonredundant~pom$AllReads), col="red")
## plot(pom$AllReads, pomNR$Nonredundant, main="RNA-Seq", xlab="All Reads", ylab="Nonredundant Reads", pch=19, xlim=c(0,3000), ylim=c(0,3000))
## abline(lm(pomNR$Nonredundant~pom$AllReads), col="red")

for (i in names(samC.70)) {
  png(filename=paste(i,".scatterplot.png", sep=""))
  plot(pom$AllReads, pomNR$Nonredundant, main=i, xlab="All Reads", ylab="Nonredundant Reads", pch=19, xlim=c(0,3000), ylim=c(0,3000), col="blue")
  points(samC.70[[i]], nrC.70[[i]], pch=19, col="black")
  abline(lm(nrC.70[[i]]~samC.70[[i]]), col="red")
  abline(lm(pomNR$Nonredundant~pom$AllReads), col="yellow")
  dev.off()
}
#dev.off()

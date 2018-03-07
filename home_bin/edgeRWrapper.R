#!/usr/bin/R

#DESeqWrapper.R by Ethan Ford

#Usage:  Rscript /NGS_scripts/DESeqWrapper.R <condition1name,condition2name> <comma separated list of conditions in same order of the columns in countstable.txt>
#Rscript edgeRWrapper.R P,C C1R5_C3,C1R5_C3_remake,C1R5_C4,C1R5_C4_remake,P1,P1_remake,P2,P2_remake
#Rscript edgeRWrapper.R P,C C,C,C,C,P,P,P,P

args <- commandArgs(TRUE)

arg1 <- unlist(strsplit(args[1], ","))
cond1 <- arg1[1]
cond2 <- arg1[2]
group <- unlist(strsplit(args[2], ","))
# cond1 <- "P"
# cond2 <- "C"
# group <- unlist(strsplit("C,C,C,C,P,P,P,P", ","))
# outname <- "test"
outname <-  gsub(",", "_", (args[1]))
library(edgeR)

# Load Data
raw.data <- read.table( file = "countstable.txt", header = TRUE )
counts <- raw.data[ , -c(1,ncol(raw.data)+1) ]
rownames( counts ) <- raw.data[ , 1 ]
# head( counts )

# Build edgeR object
## This depends on the countstable headers
#group <- c(rep("C", 4) , rep("P", 4))
cds <- DGEList( counts , group = group )

# If you need to filter out low count reads since it would be impossible to detect differential expression.
# The method used in the edgeR vignette is to keep only those genes that have at least 1 read per million in at least 3 samples.
# cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]

# normalization factors which correct for the different compositions of the samples.
cds <- calcNormFactors( cds )

pdf(paste("edgeR_",outname,".pdf",sep=""))
# To view the plot immediately
plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )

# Estimate Dispersion
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds , prior.n = 10 )

plotBCV(cds, cex=0.4)

# Mean-Variance Plot
meanVarPlot <- plotMeanVar( cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, 
    show.ave.raw.vars=FALSE, dispersion.method = "qcml", NBline = TRUE, nbins = 100, pch = 16, xlab ="Mean Expression (Log10 Scale)",
    ylab = "Variance (Log10 Scale)", main = "Mean-Variance Plot" ) 

de.aut <- exactTest( cds, pair = c( cond1 , cond2 )) # default
# de.cmn <- exactTest( cds, dispersion = "common", pair = c( "C" , "P" ) ) #single common dispersion across all genes
# de.tgw <- exactTest( cds, dispersion = "tagwise", pair = c( "C" , "P" ) ) #each gene its own dispersion
# de.poi <- exactTest( cds, dispersion = 1e-06 , pair = c( "C" , "P" ) ) #a poisson model (no dispersion)

# Store full topTags results table
resultsTbl.aut <- topTags( de.aut, n = nrow( de.aut$table ) )$table

# Names/IDs of DE genes with FDR Less than or equal to 5%
de.genes.aut <- rownames( resultsTbl.aut )[ resultsTbl.aut$FDR <= 0.05 ]

summary(de <- decideTestsDGE(de.aut, p=0.05, adjust="BH"))
detags <- rownames(cds)[as.logical(de)]
plotSmear(de.aut, de.tags=detags)
abline(h = c(-2, 2), col = "blue")

# write results
resultsTbl.aut$FC <- 2^resultsTbl.aut$logFC
resultsTbl.aut <- resultsTbl.aut[c("FC", "logFC", "logCPM", "PValue", "FDR")]

plot(resultsTbl.aut$logFC, -log10(resultsTbl.aut$PValue),
     xlim=c(-10, 10), ylim=c(0, 15), #Set limits
     main= paste("Log2 Fold Change Versus Log10 P-Value (5% FDR): ",cond1," vs ",cond2,sep=""), sub="Volcano Plot",
     col = ifelse( resultsTbl.aut$FDR < .05, "red", "black" ),
     xlab="log2 fold change", ylab="log10 p-value")#Set axis labels

#t <- is.element(rownames(resultsTbl.aut),rownames(resultsTbl.aut[resultsTbl.aut$FDR < 0.05, ]))
#points(resultsTbl.aut$logFC[t],-log10(resultsTbl.aut$PValue[t]),col='red') 

write.table(resultsTbl.aut, file=paste(outname,"_Results.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)

wh.rows.aut <- match( rownames( resultsTbl.aut ), rownames( cds$counts ) )
combResults.aut <- cbind(resultsTbl.aut, "Dispersion"=cds$tagwise.dispersion[ wh.rows.aut ], 
    "UpDown"=decideTestsDGE( de.aut, p.value = 0.05 )[ wh.rows.aut ], cds$counts[ wh.rows.aut , ] )

write.table(combResults.aut, file=paste(outname,"_ResultsDetailed.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)
write.table(cds$pseudo.alt, file=paste(outname,"_NormalizedCounts.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)

sigResultsTbl.aut <- resultsTbl.aut[ resultsTbl.aut$FDR < 0.05, ]
write.table(sigResultsTbl.aut, file=paste(outname,"_SigResults.tsv", sep=""), quote=FALSE, sep="\t", na="", col.names=NA)


library("BSgenome")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

chrset <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
uniq <- 1e-3
# Load Bam Files
cat("Loading Bams\n")
A1_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A1-HA_S5_R1_001/GSF939A-A1-HA_S5_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A2_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A2-HA_S6_R1_001/GSF939A-A2-HA_S6_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A3_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A3-HA_S7_R1_001/GSF939A-A3-HA_S7_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A4_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A4-HA_S8_R1_001/GSF939A-A4-HA_S8_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A5_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A5-HA_S9_R1_001/GSF939A-A5-HA_S9_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A6_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939A-A6-HA_S10_R1_001/GSF939A-A6-HA_S10_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A7_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A7-HA_S1_R1_001/GSF939B-A7-HA_S1_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A8_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A8-HA_S2_R1_001/GSF939B-A8-HA_S2_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A9_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A9-HA_S3_R1_001/GSF939B-A9-HA_S3_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A10_HA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A10-HA_S4_R1_001/GSF939B-A10-HA_S4_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A11_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A11-LA_S5_R1_001/GSF939B-A11-LA_S5_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A12_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A12-LA_S6_R1_001/GSF939B-A12-LA_S6_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A13_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A13-LA_S7_R1_001/GSF939B-A13-LA_S7_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A14_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A14-LA_S8_R1_001/GSF939B-A14-LA_S8_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A15_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A15-LA_S9_R1_001/GSF939B-A15-LA_S9_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A16_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF939B-A16-LA_S10_R1_001/GSF939B-A16-LA_S10_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A17_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF-939-A17-LA_S17_R1_001/GSF-939-A17-LA_S17_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A18_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF-939-A18-LA_S18_R1_001/GSF-939-A18-LA_S18_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A19_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF-939-A19-LA_S19_R1_001/GSF-939-A19-LA_S19_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)
A20_LA <- MEDIPS.createSet(file="/nfs/labs/nephew/tepper/bowtie/GSF-939-A20-LA_S20_R1_001/GSF-939-A20-LA_S20_R1_001.chr.sorted.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", shift=0, uniq=uniq, chr.select=chrset, window_size=500)

# concatenate into sets
LAset <- c(A11_LA,A12_LA,A13_LA,A14_LA,A15_LA,A16_LA,A17_LA,A18_LA,A19_LA,A20_LA)
HAset <- c(A1_HA,A2_HA,A3_HA,A4_HA,A5_HA,A6_HA,A7_HA,A8_HA,A9_HA,A10_HA)
# correlate sets
cat("Perfoming Corrleation\n")
LA.cor.matrix = MEDIPS.correlation(MSets=LAset)
HA.cor.matrix = MEDIPS.correlation(MSets=HAset)

# some quality control
# can do for all of some.
#sr_141T <- MEDIPS.saturation(file = "/nfs/labs/nephew/a2780_methylation/mbd_9-16/bowtie/merged/140911_141T.chr.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
#cr_141T <- MEDIPS.seqCoverage(file="/nfs/labs/nephew/a2780_methylation/mbd_9-16/bowtie/merged/140911_141T.chr.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", pattern = "CG", extend=250)
#er_141T <- MEDIPS.CpGenrich(file="/nfs/labs/nephew/a2780_methylation/mbd_9-16/bowtie/merged/140911_141T.chr.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19")
#MEDIPS.plotSeqCoverage(seqCoverageObj = cr_141T, type = "pie", cov.level = c(0,1,2,3,4,5))

save.image("loadedBams.Rdata")

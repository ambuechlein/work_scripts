#!/usr/bin/env R
library(methylaction2)
library("BiocParallel")

samp <- readSampleInfo("samples.csv")
chrs <- as.character(c(1:20))

fragsize <- 75
winsize <- 50
ncore <- 40 
reads <- getReads(samp=samp, chrs=chrs, fragsize=fragsize, ncore=ncore)
counts <- getCounts(samp=samp, reads=reads, chrs=chrs, winsize=winsize, ncore=ncore)

save.image("preprocess.glm.RData")

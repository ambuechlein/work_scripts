#!/bin/bash
# Compute chromosome lengths and create individual length file for each chromosome 
# The chromosome lengths are saved to a temp file - seqlens.tmp

############ seqlens.tmp ############
# 1. Chromosome
# 2. Length

## chr1	195471971
## chr2	182113224
## chr3	160039680
## chr4	156508116
## chr5	151834684
## chr6	149736546
## chr7	145441459
## chr8	129401213
## chr9	124595110
## chr10	130694993


fasta=$1
/home/mnrusimh/Scripts/General/getSequenceLengths.pl --output seqlens.tmp ${fasta}

# Each line from seqlens.tmp is saved to a chr.length file within the byChr folder
# These length files will be useful to loop through each individual chromosome

awk -F '\t' '{print > "byChr/"$1".length"}' seqlens.tmp

# Delete the tmp file

rm seqlens.tmp

#!/bin/bash
# For memory efficiency - split density files based on individual chromosomes

mkdir -p byChr
for f in *.density
do 
	grep -v '^#' ${f} | awk -F '\t' -v f=${f} '{print > "byChr/"$1"_"f}'
done

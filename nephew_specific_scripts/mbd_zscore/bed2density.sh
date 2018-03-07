#!/bin/bash
# Convert bed files to density

############ density files ############
# 1. Chromosome
# 2. Position
# 3. Coverage

## chr1	3000168	1
## chr1	3000169	1
## chr1	3000170	1
## chr1	3000171	1
## chr1	3000172	1
## chr1	3000173	1
## chr1	3000174	1
## chr1	3000175	1
## chr1	3000176	1


for f in *.bed
do
	fname=${f%.bed}
	/cluster/solexa/bin/density.pl ${f}
	mv ${f}.density ${fname}.density
done

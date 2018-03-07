#!/bin/bash
# Convert bam files to bed

############ bed format ############
# 1. Chromosome
# 2. Start
# 3. End
# 4. Read ID
# 5. MapQ
# 6. Strand

## chr1	3000167	3000203	NB500948:183:H5KH5BGX2:3:12510:6167:12899	255	-
## chr1	3000246	3000274	NB500948:183:H5KH5BGX2:4:12501:15965:11288	255	-
## chr1	3000247	3000274	NB500948:183:H5KH5BGX2:1:23209:11218:17384	255	-
## chr1	3000279	3000359	NB500948:183:H5KH5BGX2:3:11608:3771:2027	15	+
## chr1	3000374	3000454	NB500948:183:H5KH5BGX2:3:21603:7512:8814	12	-
## chr1	3000871	3000951	NB500948:183:H5KH5BGX2:4:21409:20011:15482	37	-
## chr1	3001075	3001124	NB500948:183:H5KH5BGX2:2:11312:18403:17759	30	+
## chr1	3001235	3001315	NB500948:183:H5KH5BGX2:2:13208:11579:12661	255	-
## chr1	3001726	3001806	NB500948:183:H5KH5BGX2:3:21405:8021:4064	36	-
## chr1	3003278	3003358	NB500948:183:H5KH5BGX2:4:23506:20932:19498	255	+


for f in *.bam
do
	fname=${f%.bam}
	/nfs/bio/sw/bin/bamToBed -split -i ${f} > ${fname}.bed
done

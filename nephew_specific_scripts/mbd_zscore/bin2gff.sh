# Export the regional window files to GFF format. 
# Local bins need to be matched with corresponding regional windows in order to compute Z-scores for local bins with respect to the regional windows.
# This matching is achieved using FEATnotator which requires the regional window files to be in GFF format.

############ GFF format ############
# 1. Chromosome
# 2. Source
# 3. Method
# 4. Start
# 5. End
# 6. Score
# 7. Strand
# 8. Phase
# 9. Annotation

## chr1	bin_25kx5k	bin	1	25000	.	+	.	ID=chr1_1
## chr1	bin_25kx5k	bin	5001	30000	.	+	.	ID=chr1_2
## chr1	bin_25kx5k	bin	10001	35000	.	+	.	ID=chr1_3
## chr1	bin_25kx5k	bin	15001	40000	.	+	.	ID=chr1_4
## chr1	bin_25kx5k	bin	20001	45000	.	+	.	ID=chr1_5
## chr1	bin_25kx5k	bin	25001	50000	.	+	.	ID=chr1_6
## chr1	bin_25kx5k	bin	30001	55000	.	+	.	ID=chr1_7
## chr1	bin_25kx5k	bin	35001	60000	.	+	.	ID=chr1_8
## chr1	bin_25kx5k	bin	40001	65000	.	+	.	ID=chr1_9
## chr1	bin_25kx5k	bin	45001	70000	.	+	.	ID=chr1_10

cd byChr/regional

for f in *_binned_regional.txt
do
	fname=${f%.txt}
	grep -v '^#' ${f} | awk -F '\t' '{print $1"\tbin_regional\tbin\t"$3"\t"$4"\t.\t+\t.\tID="$7}' > ${fname}.gff3
done
cd ../..

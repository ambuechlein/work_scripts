# Calculate Z difference (Z delta) between any two samples (Control and Treatment) to be compared.
# Z difference is the difference in the Z averages of the two samples

# An experimental design file listing the contrasts to be compared should be supplied as input

# The file should contain two columns separated by a tab. 
# The first column would be Control sample name and the second column would be Treatment sample name.


############ Experimental design ############
# 1. Control Sample Name
# 2. Treatment Sample Name


## Braf_ETBF_MP_Scraped	Braf_ETBF_MP_Tumor
## Braf_Sham_D_Scraped	Braf_ETBF_MP_Scraped
## Braf_Sham_MP_Scraped	Braf_ETBF_MP_Scraped
## Min_Sham_MP_Scraped	Braf_ETBF_MP_Scraped
## Min_ETBF_MP_Scraped	Braf_ETBF_MP_Scraped


############ Zdiff tables ############
# 1. Chromosome
# 2. Local Bin_Start
# 3. Local Bin_End
# 4. Local Bin Index
# 5. Regional Bin_Start
# 6. Regional Bin_End
# 7. Control Sample Z average
# 8. Treatment Sample Z average
# 9. Z difference

## chr1	1	500	chr1_1	1	25000	0	0	0
## chr1	501	1000	chr1_2	1	25000	0	0	0
## chr1	1001	1500	chr1_3	1	25000	0	0	0
## chr1	1501	2000	chr1_4	1	25000	0	0	0
## chr1	2001	2500	chr1_5	1	25000	0	0	0
## chr1	2501	3000	chr1_6	1	25000	0	0	0
## chr1	3001	3500	chr1_7	1	25000	0	0	0
## chr1	3501	4000	chr1_8	1	25000	0	0	0
## chr1	4001	4500	chr1_9	1	25000	0	0	0

experiment=$1

cd byChr

while IFS=$'\t' read -r a b
do
	control=${a}.zscore_avg
	treatment=${b}.zscore_avg
	paste ${control} ${treatment} | cut -f 1-7,14 | perl -e '
		use strict;
		use warnings;
		while (<>) {
			chomp; 
			my @cols=split /\t/; 
			if ($cols[0] eq "Chr") {
				print "$_\tZ_diff\n";
				next;
			} 
			my $zdiff = $cols[7] - $cols[6]; 
			print "$_\t$zdiff\n";
		}
	' > ${a}_x_${b}.zdiff
done < ${experiment}

cd ..


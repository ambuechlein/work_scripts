# Average Zscores across replicates for each sample.

# A sample definition file listing sample names and the corresponding replicate names should be supplied as input.
# The replicate names will be the base names of the corresponding bam files.

# The file should contain two columns separated by a tab. 
# The first column would be sample name and the second column would be a comma separated list of replicate names (no spaces)


############ Sample definition ############
# 1. Sample Name
# 2. Replicate names separated by comma

## Braf_ETBF_MP_Scraped	Braf_ETBF_MP_Scraped_5925
## Braf_ETBF_MP_Tumor	Braf_ETBF_MP_Tumor_5534_4,Braf_ETBF_MP_Tumor_5672
## Braf_Sham_D_Scraped	Braf_Sham_D_Scraped_5960,Braf_Sham_D_Scraped_6009,Braf_Sham_D_Scraped_6086
## Braf_Sham_MP_Scraped	Braf_Sham_MP_Scraped_5960,Braf_Sham_MP_Scraped_6009,Braf_Sham_MP_Scraped_6086
## Min_ETBF_MP_Scraped	Min_ETBF_MP_Scraped_5890,Min_ETBF_MP_Scraped_6073,Min_ETBF_MP_Scraped_6102
## Min_Sham_MP_Scraped	Min_Sham_MP_Scraped_5962,Min_Sham_MP_Scraped_5986,Min_Sham_MP_Scraped_6060
## Standard	Standard


############ Zaverage tables ############
# 1. Chromosome
# 2. Local Bin_Start
# 3. Local Bin_End
# 4. Local Bin Index
# 5. Regional Bin_Start
# 6. Regional Bin_End
# 7. Z average

## chr1	1	500	chr1_1	1	25000	0
## chr1	501	1000	chr1_2	1	25000	0
## chr1	1001	1500	chr1_3	1	25000	0
## chr1	1501	2000	chr1_4	1	25000	0
## chr1	2001	2500	chr1_5	1	25000	0
## chr1	2501	3000	chr1_6	1	25000	0
## chr1	3001	3500	chr1_7	1	25000	0
## chr1	3501	4000	chr1_8	1	25000	0
## chr1	4001	4500	chr1_9	1	25000	0

samples=$1

cd byChr

while IFS=$'\t' read -r a b
do
	sample=${a}
	reps=$(echo $b | tr "," "\n")
	echo -e "Chr\tStart\tEnd\tBin_Index\tRegional_Window_Start\tRegional_Window_End\t${sample}" > ${sample}.zscore_avg
	for l in *.length
	do
		c=${l%.length}
		for rep in $reps
		do
			echo ${c}_${rep}.zscores; 
		done | xargs paste | perl -e '
			use strict;
			use warnings;
			while(<>) {
				chomp; 
				my @cols=split /\t/; 
				next if ($cols[0] eq "Chr"); 
				my $ttlScore = 0; 
				my $n = (scalar @cols)/10; 
				for my $i (1..$n) {
					my $k = ($i * 10) - 1; 
					$cols[$k] = 0 if ($cols[$k] eq "NA"); 
					$ttlScore += $cols[$k];
				}
				my $avg = $ttlScore / $n; 
				print "$cols[0]\t$cols[1]\t$cols[2]\t$cols[4]\t$cols[5]\t$cols[6]\t$avg\n";
			}
		' 
	done >> ${sample}.zscore_avg 
done < ${samples}

cd ..


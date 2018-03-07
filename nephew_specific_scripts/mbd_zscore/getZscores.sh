# Compute Zscores for each local bin with respect to the corresponding overlapping regional windows.
# It is computed as the numer of standard deviations a local bin average coverage is away from an overlapping regional window average coverage.


############ Zscores tables ############
# 1. Chromosome
# 2. Local Bin_Start
# 3. Local Bin_End
# 4. Local Bin Average_Coverage
# 5. Local Bin Index
# 6. Regional Bin_Start
# 7. Regional Bin_End
# 8. Regional Window Average_Coverage
# 9. Regional Window Standard_Deviation
# 10. Z score

## chr1	1	500	0	chr1_1	1	25000	0	0	NA
## chr1	501	1000	0	chr1_2	1	25000	0	0	NA
## chr1	1001	1500	0	chr1_3	1	25000	0	0	NA
## chr1	1501	2000	0	chr1_4	1	25000	0	0	NA
## chr1	2001	2500	0	chr1_5	1	25000	0	0	NA
## chr1	2501	3000	0	chr1_6	1	25000	0	0	NA
## chr1	3001	3500	0	chr1_7	1	25000	0	0	NA
## chr1	3501	4000	0	chr1_8	1	25000	0	0	NA
## chr1	4001	4500	0	chr1_9	1	25000	0	0	NA


cd byChr

for l in *.length
do
	c=${l%.length}
	for f in ${c}_*.scores
	do
		fname=${f%.scores}
		cat ${f} | perl -e '
			use strict;
			use warnings;
			while(<>) {
				chomp; 
				my @cols=split /\t/; 
				if ($cols[0] eq "Chr") {
					print join("\t", ($cols[0], "Local_$cols[2]", "Local_$cols[3]", "Local_Bin_Coverage", "Local_Bin_Index", "Regional_Window_Start", "Regional_Window_End", "Window_Mean", "Window_StdDev", "Z_Score" )) . "\n"; 
					next;
				}
				next if ($cols[6] eq "NA");
				my $score = $cols[4]; 
				my $mean = $cols[11]; 
				my $stdev = $cols[12]; 
				my $zscore = "NA"; 
				$zscore = ($score - $mean) / $stdev unless ($stdev==0); 
				print join("\t", ($cols[0], $cols[2], $cols[3], $score, $cols[5], $cols[9], $cols[10], $mean, $stdev, $zscore)) . "\n";
			}
		' > ${fname}.zscores
	done
done

cd ..

# Match local bins with corresponding regional windows.
# It is needed in order to compute Z-scores for the local bins with respect to the regional windows.

cd byChr

# Create a workspace to process annotation information
mkdir -p annot_ws

# Run FEATnotator
for l in *.length
do
	c=${l%.length}
	for f in local/${c}_*_binned_local.txt
	do
		fname=$(basename $f)
		fname=${f%_local.txt}
		gff=regional/${fname}_regional.gff3
		/home/mnrusimh/Scripts/Annotation/annotate.pl --input ${f} --gff ${gff} --out_prefix annot_ws/${fname} --input_presorted --gff_presorted --chr_col 1 --start_col 3 --end_col 4
	done
done

# Process FEATnotator output generating tables containing matching information between local and regional windows
for l in *.length
do
	c=${l%.length}
	for f in annot_ws/${c}_*_binned.FEATnote.consolidated.txt
	do
		fname=$(basename $f)
		fname=${fname%_binned.FEATnote.consolidated.txt}
		cat ${f} | perl -e '
			use strict;
			use warnings;
			while(<>) {
				chomp;
				my @cols=split /\t/; 
				next if ($cols[0] eq "Chromosome"); 
				if ($cols[0] eq "#Chr") {
					$cols[0] = "Chr"; 
					push @cols, "Index"; 
					print join("\t", @cols) . "\n"; 
					next;
				} 
				$cols[6] =~ s/bin://g; 
				my @mappings=split /\s+/, $cols[6]; 
				foreach my $mapping (@mappings) {
					print join("\t", @cols[0..5]) . "\t$mapping\n";
				}
			}
		' > ${fname}.map
	done
done

# Generate tables listing average coverage for local bins and that along with standard deviation for regional windows


############ scores tables ############
# 1. Chromosome for Local Bin
# 2. Local Bin_ID
# 3. Local Bin_Start
# 4. Local Bin_End
# 5. Local Bin Average_Coverage
# 6. Local Bin Index
# 7. Regional Bin Index
# 8. Chromosome for Regional Window
# 9. Regional Bin_ID
# 10. Regional Bin_Start
# 11. Regional Bin_End
# 12. Regional Window Average_Coverage
# 13. Regional Window Standard_Deviation

## chr1	250	1	500	0	chr1_1	chr1_1	chr1	12500	1	25000	0	0
## chr1	750	501	1000	0	chr1_2	chr1_1	chr1	12500	1	25000	0	0
## chr1	1250	1001	1500	0	chr1_3	chr1_1	chr1	12500	1	25000	0	0
## chr1	1750	1501	2000	0	chr1_4	chr1_1	chr1	12500	1	25000	0	0
## chr1	2250	2001	2500	0	chr1_5	chr1_1	chr1	12500	1	25000	0	0
## chr1	2750	2501	3000	0	chr1_6	chr1_1	chr1	12500	1	25000	0	0
## chr1	3250	3001	3500	0	chr1_7	chr1_1	chr1	12500	1	25000	0	0
## chr1	3750	3501	4000	0	chr1_8	chr1_1	chr1	12500	1	25000	0	0
## chr1	4250	4001	4500	0	chr1_9	chr1_1	chr1	12500	1	25000	0	0


for l in *.length
do
	c=${l%.length}
	for f in ${c}_*.map
	do
		fname=${f%.map}
		/home/mnrusimh/Scripts/General/combineTwoTables.pl --ref ${f} --det regional/${fname}_binned_regional.txt --refcol 7 --detcol 7 --output ${fname}.scores
	done
done

cd ..

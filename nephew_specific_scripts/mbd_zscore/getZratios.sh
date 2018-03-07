# Calculate Z ratio from Z difference 
# Z ratio for a given local bin is the number of standard deviations the bin's Z difference is from the mean Z difference

############ Zratio tables ############
# 1. Chromosome
# 2. Local Bin_Start
# 3. Local Bin_End
# 4. Local Bin Index
# 5. Regional Bin_Start
# 6. Regional Bin_End
# 7. Control Sample Z average
# 8. Treatment Sample Z average
# 9. Z difference
# 10. Z ratio


## chr1	3531501	3532000	chr1_7064	3510001	3535000	3.98129386851136	2.96420764154394	-1.01708622696742	-4.89396056120755
## chr1	3531501	3532000	chr1_7064	3515001	3540000	3.9284857749013	2.95669941921803	-0.97178635568327	-4.6759890877825
## chr1	3531501	3532000	chr1_7064	3520001	3545000	3.94236391997561	2.99136750888482	-0.95099641109079	-4.57595315552074
## chr1	3531501	3532000	chr1_7064	3525001	3550000	3.82844691148087	2.91043017906113	-0.91801673241974	-4.4172633193417
## chr1	3531501	3532000	chr1_7064	3530001	3555000	3.80198894965918	2.89530432420769	-0.90668462545149	-4.36273609922256
## chr1	3670501	3671000	chr1_7342	3650001	3675000	-0.153388341266532	-0.223279503175188	-0.069891161908656	-0.336298517164848
## chr1	3670501	3671000	chr1_7342	3655001	3680000	-0.150661314850463	-0.208132273756738	-0.057470958906275	-0.276535655330528
## chr1	3670501	3671000	chr1_7342	3660001	3685000	-0.132963940204769	-0.186350919120976	-0.053386978916207	-0.256884581041825
## chr1	3670501	3671000	chr1_7342	3665001	3690000	-0.078424321615886	-0.17159200807473	-0.093167686458844	-0.448299240535423

cd byChr


# Calculate standard deviation of Z difference
for f in *.zdiff
do
	fname=${f%.zdiff}
	cut -f 9 ${f} | grep -v 'Z_diff' | perl -e '
		use strict;
		use warnings;
		sub mean {
			my($data) = [@_]; 
			if (not @$data) {
				die("Empty array\n"); 
			} 
			my $total = 0;
			foreach (@$data) {
				$total += $_;
			}
			my $average = $total / @$data;
			return $average;
		} 
		sub standard_deviation {
			my($data) = [@_]; 
			if(@$data == 1) { 
				return 0;
			} 
			my $average = &mean(@$data); 
			my $sqtotal = 0;
			foreach(@$data) {
				$sqtotal += ($average-$_) ** 2;
			} 
			my $std = ($sqtotal / (@$data-1)) ** 0.5; 
			return $std; 
		} 
		my @values=<>; 
		chomp @values; 
		my $stdev = &standard_deviation(@values); 
		print "'${fname}'\t$stdev\n";
	'
done > stdev.txt

# Calculate Z ratio

for f in *.zdiff
do 
	fname=${f%.zdiff}
	cat ${f} | perl -e '
		use strict;
		use warnings;
		my $stdev = 0;
		open(STDEV, "<", "stdev.txt") || die $!;
		while (<STDEV>) {
			chomp;
			my @cols=split /\t/;
			$stdev = $cols[1] if ($cols[0] eq "'${fname}'");
		} 
		while(<>) {
			chomp; 
			my @cols=split /\t/; 
			if ($cols[0] eq "Chr") {
				print "$_\tz_ratio\n";
				next;
			}
			my $zratio = $cols[8] / $stdev; 
			print "$_\t$zratio\n";
		}
	' > ${fname}.zratio
done

cd ..

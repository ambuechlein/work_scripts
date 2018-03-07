#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift;
open(my $in, $file) or die "Can't open file: $!\n";
my $h = <$in>;
chomp $h;
# 1	chr	chr1
# 2	start	12251
# 3	end	12750
# 4	width	500
# 5	pattern	Hypermethylated
# 6	all.3HH_vs_4RARA.sex_interaction.baseMean	12.4671364785581
# 7	all.3HH_vs_4RARA.sex_interaction.log2FoldChange	0.638443724245302
# 8	all.3HH_vs_4RARA.sex_interaction.lfcSE	0.268776177875975
# 9	all.3HH_vs_4RARA.sex_interaction.stat	2.37537317961232
# 10	all.3HH_vs_4RARA.sex_interaction.pvalue	0.0175312156746428
# 11	all.3HH_vs_4RARA.sex_interaction.padj	0.0301229764634748
# 12	call	intergenic
# 13	call_genes	
# 14	cpgIsland_per	0
# 15	cpgShore_per	0
# 16	cpgShelf_per	0
my @head = split(/\t/,$h);
print join("\t",@head[0,11..12,1..10,13..15]),"\n";
while(<$in>){
  chomp $_;
  my @line = split(/\t/,$_);
  print join("\t",@line[0,11..12,1..10,13..15]),"\n";
}
close $in;
exit;

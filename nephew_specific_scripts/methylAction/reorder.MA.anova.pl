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

# 0	chr	chr1
# 1	start	12101
# 2	end	12750
# 3	width	650
# 4	anodev.padj	0.00335625437502326
# 5	pattern	Hypermethylated
# 6	X3HH_vs_4RARA.p	0.00125456530335453
# 7	X3HH_over_4RARA.log2fc	0.818155268235832
# 8	call	intergenic
# 9	call_genes	
# 10	cpgIsland_per	0
# 11	cpgShore_per	0
# 12	cpgShelf_per	0

my @head = split(/\t/,$h);
print join("\t",@head[1..4,167..168,7,156,8,6,172..174]),"\n";
while(<$in>){
  chomp $_;
  my @line = split(/\t/,$_);
  print join("\t",@line[1..4,167..168,7,156,8,6,172..174]),"\n";
}
close $in;
exit;

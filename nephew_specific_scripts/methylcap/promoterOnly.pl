#!/usr/bin/perl
use strict;
use warnings;
my $file = shift;
open(my $in, $file) or die "Can't open $file: $!\n";
my $h = <$in>;
chomp $h;
my @head = split(/\t/,$h);
print "$head[0]\t",join("\t",@head[2..10]),"\n";
# Window	ENSEMBL ID Gene	ENSEMBL ID Promoter	Gene Name	Gene Type	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
# chr8:26046001-26046500	-	ENSG00000221818.8	EBF2	protein_coding promoter	17.0652158112549	-2.36184809758575	0.345998935438338	-6.82617157360202	8.72106213477073e-12	1.93070361964408e-06
# chr8:22703501-22704000	ENSG00000253125.1	-	RP11-459E5.1	processed_transcript	11.7917386602022	-2.03217975118788	0.383081971310085	-5.30481699318327	1.12786089754074e-07	0.012484517847058
# chr2:38075001-38075500	ENSG00000138061.11	-	CYP1B1	protein_coding	30.7331486212684	2.00530418533103	0.390535186262002	5.13475931458251	2.82505326559515e-07	0.0208473864050172
while(<$in>){
  chomp $_;
  my @line = split(/\t/,$_);
  if($line[2] ne '-'){
    print "$line[0]\t",join("\t",@line[2..10]),"\n";
  }
}
close $in;
exit;

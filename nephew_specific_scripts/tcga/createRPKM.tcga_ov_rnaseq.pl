#!/usr/bin/perl
use strict;
use warnings;

my $dir = shift;
my $pre = shift;
opendir(my $din, $dir) or die "Can't open $dir: $!\n";
my %rpkm;
my @order;
while (my $file = readdir($din)){
  #TCGA-04-1348-01A-01R-1565-13.hg19.gene.quantification.txt
  next unless($file =~ /(TCGA\-\d+\-\d+)\S+hg19.gene.quantification.txt$/o);
  my $id = $1;
# gene	raw_counts	median_length_normalized	RPKM
# AADACL3|126767_calculated	7	0.1297	0.0114
# AADACL4|343066_calculated	10	0.4752	0.0418
  push(@order, $id);
  open(my $in, "$dir/$file") or die "Can't open $file: $!\n";
  <$in>;
  while(<$in>){
    chomp $_;
    my @line = split(/\t/, $_);
    $rpkm{$line[0]}{$id} = $line[3];
  }
  close $in;
}
close $din;

open(my $rout, ">$pre.rpkm.tsv");
print $rout "geneID\trefID";
foreach (@order){
  print $rout "\t$_";
}
print $rout "\n";
foreach my $g (sort keys %rpkm){
  my ($x, @y) = split(/\|/o, $g);
  print $rout "$x\t", join("_", @y);
  foreach my $i (@order){
    print $rout "\t$rpkm{$g}{$i}";
  }
  print $rout "\n";
}
close $rout;
exit;

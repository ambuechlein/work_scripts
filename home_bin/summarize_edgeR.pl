#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift;
open(my $in, $file) or die "can't open file\n";
<$in>;
my $sig = 0;
my %genes;
my $up = 0;
my $down = 0;
my $nochange = 0;
print "Comaprison: $file\n";
while(<$in>){
  chomp $_;
#     0               1               2            3        4      5       6      7         8       9       10     11       12              13
# ENSEMBL ID      Gene Name       Gene Type       Chr     Strand  Start   End     FC      logFC   logCPM  PValue  FDR     Dispersion      UpDown
  my @line = split(/\t/, $_);
  next unless($line[11] <= 0.05);
  $sig++;
  $genes{$line[11]}{$line[0]} = $line[1];
  $up++ if($line[13] == 1);
  $nochange++ if($line[13] == 0);
  $down++ if($line[13] == -1);
}
close $in;
print "Total Signifincat Genes(FDR < 0.05)\t$sig\nUp\t$up\nNo Change\t$nochange\nDown\t$down\n";
my $c = 1;
foreach my $fdr (sort {$a<=>$b} keys %genes){
  foreach my $id (sort keys %{$genes{$fdr}}){
    print "$c\t$id\t$genes{$fdr}{$id}\t$fdr\n";
    $c++; 
    last if($c > 5);
  }
  last if($c > 5);
}
exit;

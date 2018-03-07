#!/usr/bin/perl
use strict;
use warnings;
my $input = shift;
my $outF = shift;
warn "Reading $input\n";
my %data;# = ();
open(my $in, $input) or die "Can't open $input: $!\n";
while(<$in>){
  chomp;
  my @cols = split /\t/;
  for my $pos ($cols[1]+1..$cols[2]) {
#    $data{$cols[0]}{$pos} = 0 unless ($data{$cols[0]}{$pos});
    $data{$cols[0]}{$pos}++;
  }
}
close $in;

open(my $out, ">$outF") or die "Cant' open $outF: $!\n";
print $out "#Ref\tPos\tCoverage\n";
foreach my $chr (sort keys %data) {
  foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}) {
    print $out "$chr\t$pos\t$data{$chr}{$pos}\n" if ($data{$chr}{$pos});
  }
}
close $out;
warn "$outF saved\n";
exit;

#!/usr/bin/env perl
use strict;
use warnings;
#my $fasta = '/research/projects/isga/tmp/junco_mappedfile.fasta';
my $fasta = shift;
my %seq;
my $header;

open(my $in, $fasta) || die "Can't open $fasta: $!\n";
foreach my $line (<$in>){
  chomp $line;
  if($line =~ /^\>/o){
    $seq{$line}++;
  }
}

foreach my $header (keys %seq){
  print "Dupe: $header\n" if $seq{$header} > 1;
}

exit;

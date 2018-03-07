#!/usr/bin/env perl
use strict;
use warnings;
my $file = shift;
open(my $in, $file) or die "can't open input\n";
open(my $out, ">$file.reduced") or die "can't open output\n";
while(<$in>){
  chomp $_;
  my @line = split(/\t/, $_);
  print $out "$line[0]\t$line[1]\t$line[2]\t$line[4]\t$line[5]\t$line[6]\t$line[13]\t$line[14]\t$line[19]\t$line[20]\t$line[21]\n";
}

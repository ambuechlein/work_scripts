#!/usr/bin/env perl
use strict;
use warnings;
my $file = shift;
my $out = $file.".new";
open(IN, $file) || die "nope\n";
open(OUT, ">$out") || die "nope2\n";
foreach my $line (<IN>){
  $line =~ s/\cM\n/\n/g;
  print OUT $line;
}
close IN;
close OUT;
exit;

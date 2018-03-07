#!/usr/bin/env perl
use strict;
use warnings;
my $f = shift;
open(my $in, $f) or die "Can't open $f: $!\n";
while(<$in>){
  chomp $_;
  my ($id, $seq) = split(/\t/,$_);
  print "$id\t";
  $seq =~ tr/atcgATCG/tagcTAGC/;
  $seq = reverse($seq);
  print "$seq\n";
}
close $in;
exit;

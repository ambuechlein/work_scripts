#!/usr/bin/env perl
use strict;
use warnings;
my $map = shift;
my $file = shift;

open(my $min, $map) or die "Can't open $map: $!\n";
#<$min>;
my %to;
my @order;
while(<$min>){
  chomp $_;
  my($f, $n) = split(/\t/, $_);
  $f =~ s/_R\d+_001//go; 
  $to{$f} = $n;
  push(@order, $n);
}
close $min;
my %print;
open(my $in, $file) or die "can't open $file: $!\n";
#<$in>;
while(<$in>){
  chomp $_;
  my @x = split(/\t/, $_);
  warn "$x[0] not in $map\n" if(not defined $to{$x[0]});
  next if(not defined $to{$x[0]});
  my $sample = $to{$x[0]};
  unshift(@x, $sample);
  $print{$sample} = join("\t", @x);
}
close $in;
print "Sample\tFilePrefix\tRawCount\tPost-TrimCount\tPercentPostTrim\n";
foreach my $s (@order){
  die "$s was not found" if(not defined $print{$s});
  print $print{$s}, "\n";
}
exit;

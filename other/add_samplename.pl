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
# GSF1536E-T32-45-4RARA-RV_S8_R1_001	4RARA	RV	male
  my($f, $n, $o, $p) = split(/\t/, $_);
  #$f =~ s/_R\d+_001//go; 
  $to{$f} = "$n\t$o\t$p";
  push(@order, $f);
}
close $min;
my %print;
open(my $in, $file) or die "can't open $file: $!\n";
#<$in>;
while(<$in>){
  chomp $_;
#  next if($_ =~ /Meth/go);
  my @x = split(/\t/, $_);
  warn "$x[0] not in $map\n" if(not defined $to{$x[0]});
  next if(not defined $to{$x[0]});
#  my $sample = $to{$x[0]};
  my $sample = $x[0];
  unshift(@x, $to{$x[0]});
  $print{$sample} = join("\t", @x);
}
close $in;
print "Sample\tFilePrefix\tRawCount\tPost-TrimCount\tPercentPostTrim\n";
foreach my $s (@order){
  #die "$s was not found" if(not defined $print{$s});
  print $print{$s}, "\n";
}
exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "dir|d=s",
                         "sampleorder|s=s",
                        );

opendir(my $din, $options{dir}) or die "Can't open $options{dir}: $!\n";
my %counts;
my @order;

my %map;
open(my $sin, $options{sampleorder}) or die "Can't open $options{sampleorder}: $!\n";
while (<$sin>) {
  chomp $_;
  my @d = split(/\t/, $_);
  $map{$d[0]} = $d[1];
  push(@order, $d[0]);
}
close $sin;

while (my $pre = readdir($din)){
# GSF974-P1-B01-PDX-1_S1_R1_001	B01-PDX-1
# GSF974-P1-B01-PDX-1_S1_R1_001/GSF974-P1-B01-PDX-1_S1_R1_001.windowCounts.tsv
  next unless(-d $pre);
  next unless( $pre =~ /^GSF/o);
  my $lab = "$pre.windowCounts.tsv";
  open(my $cin, "$options{dir}/$pre/$lab") or die "Can't open counts file for $pre/$lab: $!\n";
  warn "Opening $options{dir}/$pre/$lab\n";
  while(<$cin>){
    chomp $_;
# chr1	0	500	0
# chr1	500	1000	0
# chr1	1000	1500	0
    my ($chr, $str, $end, $c) = split(/\t/, $_);
    $str +=1;
    my $g = "$chr:$str\-$end";
    $counts{$g}{$pre} = $c;
  }
  close $cin;
}
closedir $din;

print "Window";
foreach (@order){
  print "\t$map{$_}";
}
print "\n";

foreach my $g (sort keys %counts){
  print $g;
  foreach my $p (@order){
    die "NOT DEFINED $g $p $counts{$g}{$p}\n" if(not defined $counts{$g}{$p});
    print "\t$counts{$g}{$p}";
  }
  print "\n";
}
exit;

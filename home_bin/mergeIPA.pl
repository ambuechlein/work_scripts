#!/usr/bin/env perl
use strict;
use warnings;

my $bar = shift;
my $stack = shift;

open(my $bin, $bar) or die "Can't open $bar: $!\n";
my $keep = <$bin>;
chomp $keep;
my $h1 = <$bin>;
my %merge;
my @order;
while(<$bin>){
  chomp $_;
  my @line = split(/\t/,$_);
  $merge{$line[0]} = $_;
  push(@order, $line[0]);
}
close $bin;

open(my $sin, $stack) or die "Can't open $stack: $!\n";
<$sin>;
<$sin>;
my %merged;
while(<$sin>){
  chomp $_;
  my @line = split(/\t/,$_);
  if(defined $merge{$line[0]}){
    $merge{$line[0]} .= "\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]";
  }else{
    die "$line[0] not in both\n";
  }
}
close $sin;

print $keep,"\n";
print "Ingenuity Canonical Pathways\t -log10(p-value)\tRatio\tz-score\tOverlappingMolecules (FDR < 0.05)\tDownregulated\tNo change\tUpregulated\tNo overlap with dataset\tAll Pathway Molecules\n";
foreach my $k (@order){
  print $merge{$k},"\n";
}
exit;

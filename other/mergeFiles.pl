#!/usr/bin/env perl
use strict;
use warnings;

my $f1 = shift;
my $f2 = shift;

open(my $in1, $f1) or die "Can't open $f1: $!\n";
my $h = <$in1>;
chomp $h;
my %o1;
my $c=0;
foreach my $e (split(/\t/,$h)){
  $o1{$c} = $e;
  $c++;
}
my %cnts;
my @go;
while(<$in1>){
  chomp $_;
  my @line = split(/\t/,$_);
  for(my $i=1; $i < scalar @line; $i++){
    $cnts{$line[0]}{$o1{$i}} = $line[$i];
  }
  push(@go,$line[0]);
}
close $in1;

open(my $in2, $f2) or die "Can't open $f2: $!\n";
$h = <$in2>;
chomp $h;
my %o2;
$c=0;
foreach my $e (split(/\t/,$h)){
  $o2{$c} = $e;
  $c++;
}
my $x = 0;
while(<$in2>){
  chomp $_;
  my @line = split(/\t/,$_);
  for(my $i=1; $i < scalar @line; $i++){
    $cnts{$line[0]}{$o2{$i}} += $line[$i];
  }
  die "different row order\n" if($go[$x] ne $line[0]);
  $x++;
}
close $in2;
print "GeneID";
for(my $i=1; $i < scalar keys %o1; $i++){
  print "\t$o1{$i}";
}
print "\n";

foreach my $g (@go){
  print $g;
  for(my $i=1; $i < scalar keys %o1; $i++){
    print "\t$cnts{$g}{$o1{$i}}";
  }
  print "\n";
}

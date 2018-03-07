#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "dir|d=s",
                         "map|m=s",
                        );

open(my $min, $options{map}) or die "can't open $options{map}: $!\n";
my @order;
while(<$min>){
  chomp $_;
  my @line = split(/\t/, $_);
  $line[1] =~ s/(\(|\))/_/go; $line[1] =~ s/_+/_/go;
  push(@order, $line[1]);
}
close $min;

opendir(my $din, $options{dir}) or die "Can't open $options{dir}: $!\n";
my %counts;
my %countsAll;
while (my $pre = readdir($din)){
# EZH2-T372A_S5.sorted.70bp.htseqcount.tsv
  next unless( $pre =~ /\.htseqcount\.tsv$/);
#  next if($pre =~ /deduped/o);
#  next unless($pre =~ /deduped/o);
  my $lab = $pre; $lab =~ s/\.htseqcount\.tsv//go;
  $lab =~ s/\.deduped//go;
  open(my $cin, "$pre") or die "Can't open counts file for $pre: $!\n";
  warn "Opening $pre\n";
  while(<$cin>){
    chomp $_;
    next if($_ =~ /^__/o);
    my ($g, $c) = split(/\t/, $_);
    $counts{$g}{$lab} = $c if($pre =~ /deduped/o);
    $countsAll{$g}{$lab} = $c unless($pre =~ /deduped/o);
  }
  close $cin;
}
closedir $din;
#countsAll.tsv countsDeduped.tsv
open(my $outA, ">countsAll.tsv") or die "Can't open countsAll.tsv: $!\n";
open(my $outD, ">countsDeduped.tsv") or die "Can't open countsDeduped.tsv: $!\n";

print $outA "gene";
print $outD "gene";
foreach (@order){
  print $outA "\t$_";
  print $outD "\t$_";
}
print $outA "\n";
print $outD "\n";

foreach my $g (sort keys %counts){
  print $outD $g;
  foreach my $p (@order){
    die "NOT DEFINED $g $p\n" if(not defined $counts{$g}{$p});
    print $outD "\t$counts{$g}{$p}";
  }
  print $outD "\n";
}

foreach my $g (sort keys %countsAll){
  print $outA $g;
  foreach my $p (@order){
    die "NOT DEFINED $g $p\n" if(not defined $countsAll{$g}{$p});
    print $outA "\t$countsAll{$g}{$p}";
  }
  print $outA "\n";
}
close $outA;
close $outD;
exit;

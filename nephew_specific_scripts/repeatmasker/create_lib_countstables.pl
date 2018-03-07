#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions(\%options,
                         "mapfile|m=s",
                         "repeatsdir|r=s",
                        );
my $mapfile = $options{mapfile};
my $repeatsdir = $options{repeatsdir};

$mapfile =~ /(\S+).nonredundant.map/o;
my $library = $1;
my $repeats = $repeatsdir . $library . '.nonredundant.fasta.tsv';

my %map;
open(my $in, $mapfile) or die "Can't open $mapfile: $!\n";
while(<$in>){
  chomp $_;
  my ($read, $seq) = split(/\t/, $_);
  push(@{$map{$seq}}, $read);
}
close $in;

my %rptcounts;
open(my $rpt, $repeats) or die "Can't open repeatMasker output $repeats: $!\n";
while(<$rpt>){
  next if($_ =~ /^\s*$/o);
  chomp $_;
  $_ =~ s/^\s*//o;
  my @line = split(/\s+/, $_);
  foreach my $read (@{$map{$line[4]}}){
    $rptcounts{$line[9]}++;
  }
}
close $rpt;
print "Repeat\t$library\n";
foreach my $repeat (sort keys %rptcounts){
  print "$repeat\t$rptcounts{$repeat}\n";
}
exit;

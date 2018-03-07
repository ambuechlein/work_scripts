#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $result = GetOptions(\%options,
                         "sampleorder|s=s",
                         "samlist|l=s",
                         "out|o=s",
                       );

warn "Setting Sample Order\n";
my @order;
my @altOrder;
my %map;
open(my $ofh, $options{sampleorder}) or die "Can't open $options{sampleorder}: $!\n";
while (<$ofh>) {
  chomp $_;
  my @d = split(/\t/, $_);
  push @order,$d[0];
  push @altOrder,$d[1];
  $map{$d[0]} = $d[1];
}
close $ofh;

my %genes;
open(my $sln, $options{samlist}) or die "Can't open $options{samlist}: $!\n";
while(<$sln>){
  chomp $_;
  warn "Parsing $_\n";
#  $_ =~ /(\w+).sam.gz$/o;
  $_ =~ /((\w|-)+).sam.gz$/o;
  my $pre = $1;

  die "Failed to parse prefix correctly for $_\t$pre\n" unless(defined $map{$pre});
  my $fh;
  open($fh, $_ =~ /.gz(ip)?$/ ? "zcat $_ |" : $_ =~ /.bz(ip)?2$/ ? "bzcat $_ |" : $_) || die("Open error: $_");
  while (<$fh>) {
    next unless /XT:i:1/;
    next if /XL\:Z\:tRNA/;
    if (/unique_match.(\S+)/) {
      my $match = $1;
      my @set = split /==/,$match;
      my %good = ();
      foreach my $i (@set) {
        my @d = split /~~/,$i;
        next unless $d[1];
        if ($i =~ /(\S+);/) {
          $good{$1}++;
        }
      }
      foreach my $i (keys %good) {
        $genes{$i}{$pre}++;
      }
    }
    if (/ambiguous_match.(\S+)/) {
      my $match = $1;
      my @set = split /==/,$match;
      my %good = ();
      foreach my $i (@set) {
        my @d = split /~~/,$i;
        next unless $d[1];
        if ($i =~ /(\S+);/) {
          $good{$1}++;
        }
      }
      foreach my $i (keys %good) {
        $genes{$i}{$pre}++;
      }
    }
  }
  close $fh;
}
close $sln;

open(my $fho, ">$options{out}") or die "Can't open $options{output} for writing: $!\n";
warn "Printing Output\n";
print $fho "ENSEMBLID";
foreach my $s (@order){
  print $fho "\t$map{$s}";
}
print $fho "\n";
my $t;
my $e;
foreach my $i (keys %genes) {
  $t++;
  next unless($i =~ /^ENS/);
  $e++;
  print $fho "$i";
  foreach my $s (@order){
    if(defined $genes{$i}{$s}){
      print $fho "\t$genes{$i}{$s}";
    } else {
      print $fho "\t0";
    }
  }
  print $fho "\n";
}
close $fho;

warn "Does not look like ENSEMBL ID's for Human or Mouse\nOnly ",$e/$t*100,"% of IDs were printed\n" if($e/$t < .5);

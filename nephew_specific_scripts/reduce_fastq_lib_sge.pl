#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift;
my $f = shift;
my %nr = ();

#open(my $fq, $file) or die "Can't open fastq file $file: $!\n";
my $fq;
open($fq, $file =~ /.gz(ip)?$/ ? "zcat $file |" : $file =~ /.bz(ip)?2$/ ? "bzcat $file |" : $file) || die("Open error: $file");

while(<$fq>){
  my $id = $_;
  my $seq = <$fq>;
  my $id2 = <$fq>;
  my $qual = <$fq>;
  chomp $id; chomp $seq; chomp $id2; chomp $qual;
  push(@{$nr{$seq}}, $id);
}
close $fq;

open(my $ot, ">$f.nonredundant.fasta") or die "Can't open output fasta file for $f: $!\n";
open(my $mp, ">$f.nonredundant.map") or die "Can't open output mapping file for $f: $!\n";
my $c = 1;
foreach my $seq (keys %nr){
  print $ot ">seq$c\n$seq\n";
  foreach my $r (@{$nr{$seq}}){
    $r =~ s/^\@//o;
    print $mp "$r\tseq$c\n";
  }
  $c++;
}
close $ot;
close $mp;
exit;

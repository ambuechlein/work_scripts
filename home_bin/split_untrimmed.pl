#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions(\%options,
                         "ainput|a=s",
                         "binput|b=s",
                         "outbase|o=s",
                        );

#open(my $ain, $options{ainput}) or die "can't open first fastq file $options{ainput}: $!\n";
my $ain;
open($ain, $options{ainput} =~ /.gz(ip)?$/ ? "zcat $options{ainput} |" : $options{ainput} =~ /.bz(ip)?2$/ ? "bzcat $options{ainput} |" : $options{ainput}) || die("Open error: $options{ainput}");
my %reads;
print STDOUT "Parsing file A\n";
while(<$ain>){
  chomp $_;
  my $seq = <$ain>;
  my $q1 = <$ain>; 
  my $q2 = <$ain>;
  chomp $seq;
  $reads{$_} = $seq;
}
close $ain;
print STDOUT "Total file A reads:\t", scalar keys %reads, "\n";

open(my $bin, $options{binput}) or die "can't open second fastq file $options{binput}: $!\n";
open(my $uout, ">$options{outbase}.untrimmed.fastq") or die "can't open output file untrimmed: $!\n";
open(my $tout, ">$options{outbase}.trimmed.fastq") or die "can't open output file trimmed: $!\n";
print STDOUT "Parsing file B and creating $options{outbase}.untrimmed.fastq and $options{outbase}.trimmed.fastq files\n";
while(<$bin>){
  chomp $_;
  my $seq = <$bin>;
  my $q1 = <$bin>;
  my $q2 = <$bin>;
  chomp $seq;
  die "Read not previously defined: $_\n" if(not defined $reads{$_});
  if ($seq eq $reads{$_}){
# write to untrimmed
    print $uout "$_\n$seq\n".$q1.$q2;
  } else {
# write to trimmed
    print $tout "$_\n$seq\n".$q1.$q2;
  }
}
close $ain;


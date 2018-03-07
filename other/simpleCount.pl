#!/usr/bin/env perl
use strict; 
use warnings;
use Getopt::Long;
use Data::Dumper;
my %options;
my $results = GetOptions( \%options,
                          "afastq|a=s",
                        );
countFastq($options{afastq});
exit;

sub countFastq{
  my $input = shift;
  my $in;
  open($in, $input =~ /.gz(ip)?$/ ? "zcat $input |" : $input =~ /.bz(ip)?2$/ ? "bzcat $input |" : $input) || die("Open error: $input");
  my $count = 0;
  while (<$in>) {
      my $line0= $_;
      my $line1=<$in>;
      my $line2=<$in>;
      my $line3=<$in>;
      $count++;
  }
  my $prefix = $input;
  $prefix =~ s/\.fastq(\.gz)*//go;
  print "$prefix\t$count\n";
  close $in;
}

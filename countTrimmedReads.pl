#!/usr/bin/env perl
use strict; 
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions( \%options,
                          "afastq|a=s",
                          "bfastq|b=s",
                          "output|o=s",
                          "prefix|p=s",
                        );
my $c1 = countFastq($options{afastq});
my $c2 = countFastq($options{bfastq});
my $survived = $c2/$c1 * 100;

$c1 =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$c2 =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;

my $prefix = $options{prefix};
$prefix =~ s/\.gz$//go;
$prefix =~ s/\.fastq$//go;
open(my $out, ">$options{output}") or die "Can't open output file: $!\n";
print $out "$prefix\t$c1\t$c2\t$survived\n";
close $out;
exit;

sub countFastq{
  my $input = shift;
  open(INPUT, $input =~ /.gz(ip)?$/ ? "zcat $input |" : $input =~ /.bz(ip)?2$/ ? "bzcat $input |" : $input) || die("Open error: $input");
  my $count = 0;
  while (<INPUT>) {
      my $line0= $_;
      my $line1=<INPUT>;
      my $line2=<INPUT>;
      my $line3=<INPUT>;
      $count++;
  }
  close INPUT;
  
  print "$input\t$count\n";
  return $count;
}

#!/usr/bin/env perl
use strict; 
use warnings;
use Getopt::Long;
my $input;
GetOptions(
          'fastq=s'  => \$input,
          );
open(INPUT, $input =~ /.gz(ip)?$/ ? "zcat $input |" : $input =~ /.bz(ip)?2$/ ? "bzcat $input |" : $input) || die("Open error: $input");
my $count = 0;
while (<INPUT>) {
    my $line0= $_;
    my $line1=<INPUT>;
    my $line2=<INPUT>;
    my $line3=<INPUT>;
    $count+=length($line1);
}
close INPUT;

print "$input\t$count\n";
exit;

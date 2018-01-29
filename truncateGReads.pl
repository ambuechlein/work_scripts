#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $result = GetOptions(\%options,
                        "input|i=s",
                        "output|o=s",
                        "length|l=s",
                       );
my $input = $options{input};
my $in;
open($in, $input =~ /.gz(ip)?$/ ? "zcat $input |" : $input =~ /.bz(ip)?2$/ ? "bzcat $input |" : $input) || die("Open error: $input");
open(my $out, ">$options{output}") or die "Can't open output fastq $options{output}: $!\n";
while (<$in>) {
    my $line0= $_;
    my $line1=<$in>;
    my $line2=<$in>;
    my $line3=<$in>;
    chomp $line0;
    chomp $line1;
    chomp $line2;
    chomp $line3;
    if($line1 =~ /GGG$/o){
      $line1 =~ s/G+$//gio;
      my $sl = length($line1);
      $line3 = substr($line3, 0, $sl);

    }
    next if(length($line1) < $options{length});
    print $out "$line0\n$line1\n$line2\n$line3\n";
}
close $in;
close $out;
exit;


#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $result = GetOptions(\%options,
                        "in|i=s",
                        "out|o=s",
                       );
my $in;
open($in, "samtools view -h $options{in} |");

#open($in, $options{in} =~ /.gz(ip)?$/ ? "zcat $options{in} |" : $options{in} =~ /.bz(ip)?2$/ ? "bzcat $options{in} |" : $options{in}) || die("Open error: $options{in}: $!\n");
open(my $out, ">$options{out}") or die "Can't open output $options{out}: $!\n";
while(<$in>){
  chomp $_;
  if($_ =~ /^\@/){
    print $out "$_\n";
    next;
  }
  my @x = split(/\t/, $_);
  print $out "$_\n" if(length($x[9]) > 45);
}
close $in;
exit;

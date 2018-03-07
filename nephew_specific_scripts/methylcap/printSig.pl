#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $result = GetOptions(\%options,
                        "sig|s=s",
                        "table|t=s",
                       );

open(my $tin, $options{table}) or die "Can't open input $options{table}: $!\n";
my $h = <$tin>;
print $h;
while(<$tin>){
  chomp $_;
  my @line = split(/\t/, $_);
  my $fdr = pop @line;
  next if($fdr eq 'NA');
  print "$_\n" if($fdr < $options{sig});
}
close $tin;
exit;

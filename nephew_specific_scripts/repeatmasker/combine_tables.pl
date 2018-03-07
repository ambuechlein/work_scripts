#!/usr/bin/env perl
use strict;
use warnings;

my %countstable;
my @pools;
foreach my $file ( @ARGV ){
  open(my $in, $file) or die "Can't open file $file: $!\n";
  my $header = <$in>;
  chomp $header;
  my ($r, $pool) = split(/\t/, $header);
  push(@pools, $pool);
  while(<$in>){
    chomp $_;
    my($repeat, $count) = split(/\t/, $_);
    $countstable{$repeat}{$pool} = $count;
  }
  close $in;
}

print "Repeat";
foreach (sort @pools) {print "\t$_";}
foreach my $repeat (sort keys %countstable){
  print "\n$repeat";
  foreach my $pool (sort @pools){
    if (defined $countstable{$repeat}{$pool}){
      print "\t$countstable{$repeat}{$pool}";
    } else {
      print "\t0";
    }
  }
}
exit;

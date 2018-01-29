#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my $input;

GetOptions(
        "i|input:s"    => \$input,
);

open(my $in, $input) or die "Can't open $input: $!\n";
#my $name_col=0;
#my $score_col=13;
my %hits;
while(<$in>){
  s/\r?\n//;
  my @F=split /\t/, $_;
  my ($id, $hit_id, $bitscore) = @F[0, 5, 13];
  if(not defined $hits{$id}{$bitscore}){
    $hits{$id}{$bitscore} = $hit_id;
  } else {
    $hits{$id}{$bitscore} .= "\n$hit_id";
  }
}
close $in;

foreach my $id (keys %hits){
  my $c = 0;
  foreach my $bitscore (sort {$b <=> $a} keys %{$hits{$id}}){
    print "$hits{$id}{$bitscore}\n";
    $c++;
    last if($c > 10);
  }
}
exit

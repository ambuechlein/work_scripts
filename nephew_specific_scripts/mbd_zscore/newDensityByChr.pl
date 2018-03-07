#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
my %options;
my $result = GetOptions(\%options,
                        "input|i=s",
                        "outF|o=s",
                        "split|s=s",
                       );

my $input = $options{input};
my $outF = $options{outF};
my $split = $options{split};
$split = 0 unless(defined $split);
my ($outB, $dir) = fileparse($outF);
system("mkdir -p $dir/byChr") if($split);
warn "Reading $input\n";
my %data = ();
open(my $in, $input) or die "Can't open $input: $!\n";
open(my $out, ">$outF") or die "Cant' open $outF: $!\n";
print $out "#Ref\tPos\tCoverage\n";
my $chr;
my %seen;
while(<$in>){
  chomp;
  my @cols = split /\t/;
  $chr = $cols[0] if(not defined $chr);
  if($chr eq $cols[0]){
    for my $pos ($cols[1]+1..$cols[2]) {
      $data{$pos}++;
    }
  } else {
    die "Bed not ordered, please order bed by position\n" if(defined $seen{$chr} and $seen{$chr} > 1);
    warn "Writing $chr\n";
    my $chrF = "$dir/byChr/${chr}_$outB";
    open(my $outC, ">$chrF") or die "Can't open $chrF: $!\n" if($split);
    foreach my $pos (sort {$a <=> $b} keys %data) {
      print $out "$chr\t$pos\t$data{$pos}\n" if ($data{$pos});
      print $outC "$chr\t$pos\t$data{$pos}\n" if($split);
    }
    close $outC if($split);
    $seen{$chr}++;
    %data = ();
    $chr = $cols[0];
    warn "Reading $chr\n";
    for my $pos ($cols[1]+1..$cols[2]) {
      $data{$pos}++;
    } 
  }
}
close $in;
exit;

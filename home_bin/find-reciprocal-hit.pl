#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my ($a_file, $b_file, $outF);

GetOptions(
        "a|a_file:s"    => \$a_file,
        "b|b_file:s"   => \$b_file,
        "o|out:s" => \$outF,
);

my %best_a_b = &find_bbh($a_file);
my %best_b_a = &find_bbh($b_file);
#
my $col1=1;
my $col2=0;
my %line1;
my @aba;
while ((my $key, my $value) = each(%best_a_b)) {
  $value =~ s/\r?\n//;
  my @F=split /\t/, $value;

  $line1{$F[$col1]} .= "$value\n";
}

while ((my $key, my $value) = each(%best_b_a)) {
  $value =~ s/\r?\n//;
  my @F=split /\t/, $value; 
  if (my $x = $line1{$F[$col2]}) {
    $x =~ s/\n/\t$value\n/g; 
    push @aba, $x;
  }
}
#
my $colm=0; 
my $coln=13;
my $count=0;
my @aba_recip;
my $totC=0;
foreach(@aba) {
  s/\r?\n//; 
  my @F=split /\t/, $_; 
  if ($F[$colm] eq $F[$coln]) {
    push @aba_recip, "$_\n"; 
    $count++;
  }
  $totC++;
} 
warn "\nChose $count lines out of $totC where column $colm had same text as column $coln\n";
#
my @cols=(0, 1, 11);
$totC=0;
open(my $out, ">$outF") or die "Can't open $outF:$!\n";
foreach(@aba_recip) {
  s/\r?\n//; 
  my @F=split /\t/, $_; 
  print $out join("\t", @F[@cols]), "\n";
  $totC++;
} 
warn "\nJoined columns ", join(", ", @cols), " for $totC lines\n";
#
exit;

sub find_bbh{
  my $file = shift;
  open(my $a_b, $file) or die "Can't open $file: $!\n";
  my $name_col=0;
  my $score_col=11;
  my (%max, %best, @names);
  while(<$a_b>){
    s/\r?\n//;
    my @F=split /\t/, $_;
    my ($n, $s) = @F[$name_col, $score_col];
    push @names, $n if (! exists($max{$n}));
    if (! exists($max{$n}) || $s > $max{$n}) {
      $max{$n} = $s;
      $best{$n} = ();
    }
    if ($s == $max{$n}) {
      $best{$n} .= "$_\n";
    } 
  }
  close $a_b;
  return %best;
}

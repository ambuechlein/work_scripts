#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my ($est_vs_orthodb, $orthodb_vs_est, $rbbh, $helpFlag);

GetOptions(
        "h|helpFlag:s"  => \$helpFlag,
        "e|est_vs_orthodb:s"    => \$est_vs_orthodb,
        "o|orthodb_vs_est:s"   => \$orthodb_vs_est,
        "r|rbbh:s" => \$rbbh,
);

my %best_e_o = &find_bbh($est_vs_orthodb);
#print Dumper(%best_e_o); exit;
my %best_o_e = &find_bbh($orthodb_vs_est);
#print Dumper(%best_o_e); exit;

#
my $col1=5;
my $col2=0;
my %line1;
my @aba;
while ((my $key, my $value) = each(%best_e_o)) {
  $value =~ s/\r?\n//;
  my @F=split /\t/, $value;
  $line1{$F[$col1]} .= "$value\n";
#print "$F[$col1]\n$line1{$F[$col1]}\n$value\n"; exit;
}

while ((my $key, my $value) = each(%best_o_e)) {
  $value =~ s/\r?\n//;
  my @F=split /\t/, $value; 
  if (my $x = $line1{$F[$col2]}) {
    $x =~ s/\n/\t$value\n/g; 
    push @aba, $x;
#print "$x\n", Dumper(@aba); exit;
  }
}
#
my $colm=0; 
my $coln=26;
my $count=0;
my @aba_recip;
foreach(@aba) {
  s/\r?\n//; 
  my @F=split /\t/, $_;
# print "$F[$colm]\t$F[$coln]\n", Dumper(@F), "\n"; exit;
  if ($F[$colm] eq $F[$coln]) {
    push @aba_recip, "$_\n"; 
    $count++;
  }
} 
warn "\nChose $count lines out of ", scalar @aba, "  where column $colm had same text as column $coln\n\n";
#
my @cols=(0, 5, 13);
open(my $out, ">$rbbh") or die "Can't open $rbbh:$!\n";
print $out "$est_vs_orthodb\t$orthodb_vs_est\tbitscore\n";
foreach(@aba_recip) {
  s/\r?\n//; 
  my @F=split /\t/, $_; 
  print $out join("\t", @F[@cols]), "\n";
} 
warn "\nJoined columns ", join(", ", @cols), " for ", scalar @aba_recip, " lines\n\n";
#
exit;

sub find_bbh{
  my $file = shift;
  open(my $a_b, $file) or die "Can't open $file: $!\n";
  my $name_col=0;
  my $score_col=13;
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

## open(my $b_a, $orthodb_vs_est) or die "Can't open $orthodb_vs_est: $!\n";
## my (%max, %best);
## while(<$b_a>){
##   /\r?\n//;
##   my @F=split /\t/, $_;
##   my ($n, $s) = @F[$name_col, $score_col];
##   push @names, $n if (! exists($max{$n}));
##   if (! exists($max{$n}) || $s > $max{$n}) {
##     $max{$n} = $s;
##     $best{$n} = ();
##   }
##   if ($s == $max{$n}) {
##     $best{$n} .= "$_\n";
##   }
##   for $n (@names) {
##     print $best{$n}
##   }
## }
## close $b_a;
#################################
## my $col1=1;
## my $col2=0;
## while (<F1>) {
##   s/\r?\n//; 
##   @F=split /\t/, $_; 
##   $line1{$F[$col1]} .= "$_\n";
## }
## 
## while (<F2>) {
##   s/\r?\n//;@F=split /\t/, $_; 
##   if ($x = $line1{$F[$col2]}) {
##     $x =~ s/\n/\t$_\n/g; 
##     print $x
##   }
## }
#################################
## $colm=0; 
## $coln=13; 
## $count=0;
## while(<>) {
##   s/\r?\n//; 
##   @F=split /\t/, $_; 
##   if ($F[$colm] eq $F[$coln]) {
##     print "$_\n"; $count++;
##   }
## } 
## warn "\nChose $count lines out of $. where column $colm had same text as column $coln\n\n";
#################################
## my @cols=(0, 1, 11); 
## while(<>) {
##   s/\r?\n//; 
##   @F=split /\t/, $_; 
##   print join("\t", @F[@cols]), "\n";
##   }
## warn "\nJoined columns ", join(", ", @cols), " for $. lines\n\n";
#################################

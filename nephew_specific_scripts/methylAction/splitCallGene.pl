#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;
open(my $in, $file) or die "Can't open $file: $!\n";
my $h = <$in>;
chomp $h;
my @head = split(/\t/,$h);
# 31	call
my $c = 0;
my $callCol;
foreach my $col (@head){
  if($col eq "call_genes"){
    $callCol = $c;
  }
  $c++;
}

print "$h\n";

while(<$in>){
  chomp $_;
  my @line = split(/\t/,$_);
  if($line[$callCol] =~ /\,/go){
    my @genes = split(/\,/,$line[$callCol]);
    foreach my $g (@genes){
      $line[$callCol] = $g;
      print join("\t",@line), "\n";
    }
  } else {
    print join("\t",@line), "\n";
  }
}
exit;

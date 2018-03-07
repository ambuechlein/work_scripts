#!/usr/bin/env perl
use strict;
use warnings;
my $file = shift;
open(my $in, $file) or die "Can't open $file: $!\n";
open(my $ot, ">$file.tmp") or die "Can't open $file.tmp: $!\n";
while(<$in>){
  chomp $_;
  my @line = split(/\t/,$_);
  print $ot "$_\n" if($line[3] =~ /\/2$/o);
  next if($line[3] =~ /\/2$/o);
  if($line[5] eq '-'){
    $line[5] = '+';
  }elsif($line[5] eq '+'){
    $line[5] = '-';
  }else{
    die "Bad BAM file: $_\n";
  }
#  print $ot join(@line,"\t"),"\n";
  print $ot join("\t",@line),"\n";
}
close $in;
close $ot;
rename("$file.tmp",$file);
exit;

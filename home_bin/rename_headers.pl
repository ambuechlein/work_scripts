#!/usr/bin/env perl
use strict;
use warnings;
my $file = shift;
open (my $in, $file) or die "can't open in file: $!\n";
open (my $out, ">$file.new") or die "can't open out: $!\n";
my %seen;
my $c = 0;
foreach my $line (<$in>){
  chomp $line;
  if($line =~ /^\>(\S+)(.*)/){
    $c++;
    if($seen{$1}){
      my $tmp = $1;
      until (not $seen{$tmp}){
        $tmp .= "A";
      }
#      my $header = '>'.$1."A".$2;
      my $header = '>'.$tmp.$2;

      print $out "$header\n";
      print "$line\n$header\n\n";
      $seen{$tmp} = 1;

    }else{
      $seen{$1} = 1;
      print $out "$line\n"
    }
  }else{
    print $out "$line\n";
  }

}
print $c, "\n";
close $in;
close $out;
exit;

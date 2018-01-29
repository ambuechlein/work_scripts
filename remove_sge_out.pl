#!/usr/bin/env perl
use strict;
use warnings;

foreach my $file (@ARGV){
  my $line = `tail -n 1 $file`;
  chomp $line;
  if($line eq 'Finished'){
    $file =~ /(\S+)\.sh\.o(\d+)/o;
    my $file2 = $1.'.sh.e'.$2;
    my $file3 = $1.'.sh.pe'.$2;
    my $file4 = $1.'.sh.po'.$2;
    if(-e $file3 and -e $file4){
      system("rm $file $file2 $file3 $file4");
    } else {
      system("rm $file $file2");
    }
  }
}
exit;

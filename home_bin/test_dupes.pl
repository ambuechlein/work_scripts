#!/usr/bin/env perl
use strict;
use warnings;
#my $fasta = '/research/projects/isga/tmp/junco_mappedfile.fasta';
my $fasta = '/nfs/bio/db/OrthoDB/OrthoDB4/fasta/OrthoDB';
my %seq;
my ($header, $est);
#my $test1 = '>Chr8:18080029-18080086 +';
#my $test2 = '>Chr14:8454572-8454648 +';
open(my $in, $fasta) || die "Can't open $fasta: $!\n";
foreach my $line (<$in>){
  chomp $line;
  if($line =~ /^\>/o){
    if(defined $est and $est ne ''){
      $seq{$header}{$est}++;
    }
    $header = $line;
    $est = '';
  }else{
    $est .= $line;
  }
}
$seq{$header}{$est}++;
foreach my $header (keys %seq){
  my $uniq = keys %{$seq{$header}};
  print "$header\n" if $uniq > 1;
}
my $tot = keys %seq;
print "Total Sequences: $tot\n\n";

exit;

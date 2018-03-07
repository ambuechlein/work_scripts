#!/usr/bin/env perl
use strict;
use warnings;

#my $fasta = '/research/projects/isga/tmp/junco_mappedfile.fasta';
#my $new = '/research/projects/isga/tmp/junco_mappedfile.fasta.new';
#my $new = '/home/abuechle/junco_mappedfile.fasta.new';
my $fasta = '/research/projects/isga/prod/project/output_repository/OrfPredictor/146155620_default/i1/g1/concatenate_files.OrfPredictor.res';
my $new = '/research/projects/isga/prod/project/output_repository/OrfPredictor/146155620_default/i1/g1/concatenate_files.OrfPredictor.res.new';

my %seq;
my ($header, $est);

open(my $in, $fasta) || die "Can't open $fasta: $!\n";
open(my $out, ">$new") || die "Can't open $new: $!\n";

foreach my $line (<$in>){
  chomp $line;
  if($line =~ /^\>/o){
    if(defined $est and $est ne ''){
      $seq{$header}{$est}++;
      print $out "$header\n$est\n" if($seq{$header}{$est} == 1);
    }
    $header = $line;
    $est = '';
  }else{
    $est .= $line;
  }
}

my $tot = keys %seq;
print "Total Sequences: $tot\n";
exit;

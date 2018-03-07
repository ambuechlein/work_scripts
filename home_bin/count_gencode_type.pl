#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "gencode_file|g=s",
                        );
my $gencode_file = $options{gencode_file};
open(my $in, "<$gencode_file") or die "Can't open $gencode_file.\n";

my %gene_types;
while(<$in>){
  next if(/^#/); #ignore header
  chomp;
  my ($chr, $source, $type, $start, $end, $score, 
  $strand, $phase, $attributes) = split("\t");
  next unless($type eq 'gene');
  my @add_attributes = split(";", $attributes);
  foreach my $attr ( @add_attributes ) {
# gene_id "ENSG00000228630.1"; transcript_id "ENSG00000228630.1"; gene_type "antisense"; gene_status "KNOWN"; gene_name "HOTAIR"; transcript_type "antisense"; transcript_status "KNOWN"; transcript_name "HOTAIR"; level 2; havana_gene "OTTHUMG00000152934.1";

#     next unless $attr =~ /^\s*(.+)\s(.+)$/;
     next unless $attr =~ /^\s*(gene_type)\s(.+)$/;
#     next unless $attr =~ /^\s*(gene_biotype)\s(.+)$/;
     my $c_type  = $1;
     my $c_value = $2;
     $c_value =~ s/\"//go;
     $gene_types{$c_value}++;
  }
}
close $in;

foreach my $type (sort {$gene_types{$b}<=>$gene_types{$a}} keys %gene_types){
  print "$type\t$gene_types{$type}\n";
}
exit;

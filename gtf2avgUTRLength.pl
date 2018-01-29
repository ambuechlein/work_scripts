#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util 'sum';
use Math::Round;
my %options;
my $results = GetOptions(\%options,
                         "gencode_file|g=s",
                         "output|o=s",
                        );
my $gencode_file = $options{gencode_file};
open(my $in, "<$gencode_file") or die "Can't open $gencode_file.\n";
my %lengths;
while(<$in>){
  next if(/^\#/); #ignore header
  chomp;
  my ($chr, $source, $type, $start, $end, $score, 
    $strand, $phase, $attributes) = split("\t");
  next unless($type eq 'UTR');
  #store nine columns in hash
  my %fields = (
    chr        => $chr,
    source     => $source,
    type       => $type,
    start      => $start,
    end        => $end,
    score      => $score,
    strand     => $strand,
    phase      => $phase,
    attributes => $attributes,
  );
  my @add_attributes = split(";", $attributes);
  # store ids and additional information in second hash
  foreach my $attr ( @add_attributes ) {
     next unless $attr =~ /^\s*(.+)\s(.+)$/;
     my $c_type  = $1;
     my $c_value = $2;
     $c_value =~ s/\"//go;
     if($c_type  && $c_value){
       $fields{$c_type} = $c_value;
     }
  }
  die "UTR does not have a gene_id: $attributes\n" if (not defined $fields{gene_id});
  die "UTR does not have a transcript_id: $attributes\n" if (not defined $fields{transcript_id});
  my $l = $fields{end} - $fields{start};
#  $lengths{$fields{gene_id}}{$fields{transcript_id}} += $l;
  $lengths{$fields{transcript_id}} += $l;

}
close $in;
#print Dumper(%lengths);
open(my $out, ">$options{output}") or die "Can't open output: $!\n";
print $out "\#GeneID\tAvgTranscriptLength\n";
foreach my $gene (sort keys %lengths){
#  my $val = sum values %{$lengths{$gene}};
#  my $kl = scalar keys %{$lengths{$gene}};
#  print "$val\n";
  warn "Length: ", $lengths{$gene}, "\n";
#  my $avg = $val/$kl;
#  my $ravg = round($avg);
#  my $test = join("\t", values %{$lengths{$gene}});
#  print $out "$gene\t$ravg\n";
  print $out "$gene\t $lengths{$gene}\n";

}
exit;

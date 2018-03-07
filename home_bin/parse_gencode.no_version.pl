#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "gencode_file|g=s",
                         "table|t=s",
                        );
my $gencode_file = $options{gencode_file};
my $table = $options{table};
open(my $in, "<$gencode_file") or die "Can't open $gencode_file.\n";

my %all_genes;

while(<$in>){
  next if(/^\#/); #ignore header
  chomp;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score, 
  $strand, $phase, $attributes) = split("\t");
  next unless($type eq 'gene');
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
       if(!exists($attribs{$c_type})){
         $attribs{$c_type} = [];
       }
       push(@{ $attribs{$c_type} }, $c_value);
       $fields{$c_type} = $c_value;
     }
  }

  #work with the information from the two hashes...
  #eg. store them in a hash of arrays by gene_id:
  if(not defined $attribs{'gene_id'}){
#    warn Dumper(@add_attributes); 
    next;
  }
#  warn $attribs{'gene_id'}->[0], "\n";
# 1       pseudogene      gene    3054233 3054733 .       +       .       gene_id "ENSMUSG00000090025"; gene_name "Gm16088"; gene_source "havana"; gene_biotype "pseudogene";
  $attribs{'gene_id'}->[0] =~ s/(ENSMUSG\d+)\.\d+/$1/o;
#  die $attribs{'gene_id'}->[0], "\n";
  if(!exists($all_genes{$attribs{'gene_id'}->[0]})){
    $all_genes{$attribs{'gene_id'}->[0]} = [];
  }
  push(@{ $all_genes{$attribs{'gene_id'}->[0]} }, \%fields);
}
close $in;
#print ${${$all_genes{'ENSG00000223972.4'}}[0]}{"type"}, "\n";
open(my $tb, "$table") or die "Can't open table file: $!\n";
my $header = <$tb>;
my @h = split(/\t/, $header);
my $hid = shift @h;
my $print_h = join("\t", "ENSEMBL ID","Gene Name","Gene Type","Chr","Strand","Start","End",@h);
print "$print_h";

while(<$tb>){
  my @line = split(/\t/, $_);
  my $id = shift @line;
  $id =~ s/(ENSG\d+)\.\d+/$1/o;
  die "$id not defined\n" if(not defined $all_genes{$id});
  my $gene_name = ${${$all_genes{$id}}[0]}{"gene_name"};
  my $chr = ${${$all_genes{$id}}[0]}{"chr"};
  my $start = ${${$all_genes{$id}}[0]}{"start"};
  my $end = ${${$all_genes{$id}}[0]}{"end"};
  my $strand = ${${$all_genes{$id}}[0]}{"strand"};
#  my $type = ${${$all_genes{$id}}[0]}{"gene_type"};
  my $type = ${${$all_genes{$id}}[0]}{"gene_biotype"};
  my $print = join("\t", $id,$gene_name,$type,$chr,$strand,$start,$end,@line);
  print $print;

}

exit;


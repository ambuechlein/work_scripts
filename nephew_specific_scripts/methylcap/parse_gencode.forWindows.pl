#!/usr/bin/perl
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
  next if(/^##/); #ignore header
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
  if(!exists($all_genes{$attribs{'gene_id'}->[0]})){
    $all_genes{$attribs{'gene_id'}->[0]} = [];
  }
  push(@{ $all_genes{$attribs{'gene_id'}->[0]} }, \%fields);
}
close $in;
#print ${${$all_genes{'ENSG00000223972.4'}}[0]}{"type"}, "\n";
open(my $tb, "$table") or die "Can't open table file: $!\n";
my $header = <$tb>;
chomp $header;
my @h = split(/\t/, $header);
#my $hid = shift @h;
#my $print_h = join("\t", "ENSEMBL ID","Gene Name","Gene Type","Chr","Strand","Start","End",@h);
shift @h; shift @h; shift @h;
my $print_h = join("\t", "Window","ENSEMBL ID Gene","ENSEMBL ID Promoter","Gene Name","Gene Type",@h);
print "$print_h\n";

while(<$tb>){
  chomp $_;
  my @line = split(/\t/, $_);
  my $wid = shift @line;
  my $gid = shift @line;
  my $pid = shift @line;
  my $gene_name = "-";
  my $type = "-";
  if($gid ne "-"){
    $gene_name = ${${$all_genes{$gid}}[0]}{"gene_name"};
#  my $chr = ${${$all_genes{$id}}[0]}{"chr"};
#  my $start = ${${$all_genes{$id}}[0]}{"start"};
#  my $end = ${${$all_genes{$id}}[0]}{"end"};
#  my $strand = ${${$all_genes{$id}}[0]}{"strand"};
    $type = ${${$all_genes{$gid}}[0]}{"gene_type"};
  }
  if($pid ne "-"){
    $gene_name = ${${$all_genes{$pid}}[0]}{"gene_name"};
    $type = ${${$all_genes{$pid}}[0]}{"gene_type"} . " promoter";
  }

  my $print = join("\t", $wid,$gid,$pid,$gene_name,$type,@line);
  print $print, "\n";

}

exit;


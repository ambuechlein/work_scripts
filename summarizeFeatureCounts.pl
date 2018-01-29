#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "gencode_file|g=s",
                         "dir|d=s",
                        );
my $gencode_file = $options{gencode_file};
my $table = $options{table};
warn "Parsing $gencode_file\n";
open(my $in, "<$gencode_file") or die "Can't open $gencode_file.\n";
my %all_genes;

while(<$in>){
  next if(/^##/); #ignore header
  chomp;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score, 
  $strand, $phase, $attributes) = split("\t");
  next unless($type eq 'gene' or ($type eq 'exon' and $source eq 'hg38_rmsk') or ($chr !~ /^chr/ and $type eq 'exon') );
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

  if(not defined $attribs{'gene_id'}){
    next;
  }
  if(!exists($all_genes{$attribs{'gene_id'}->[0]})){
    $all_genes{$attribs{'gene_id'}->[0]} = [];
  }
  push(@{ $all_genes{$attribs{'gene_id'}->[0]} }, \%fields);
}
close $in;

opendir(my $din, $options{dir}) or die "Can't open $options{dir}: $!\n";
my @samp;
my %hash;
while(my $file = readdir($din)){
  next unless($file =~ /.featureCountsaccepted_hits.sorted.bam$/o);
  warn "Parsing $file\n";
  my $pre = $file;
  $pre =~ s/\.featureCountsaccepted_hits.bam//go;
  $pre =~ s/\.accepted_hits\S*\.bam//go;
  $pre =~ s/_thout//go;
  $pre =~ s/\.fastq.gz//go;
  push(@samp,$pre);
  open(my $in, $file) or die "Can't open $file: $!\n";
  my %seen = ();
  while(<$in>){
    chomp $_;
#   NB500948:227:H3WFCBGX3:4:23609:2002:8792	Assigned	1	MER11C
    my ($read, $assign, $num, $id) = split(/\t/,$_);
    my $type = 'TotalMappingsUnassigned';
if($num <= 1){
    $type = ${${$all_genes{$id}}[0]}{"gene_type"} if($id ne '*' and $id ne 'NA');
    die "$_\n" if(not defined $type);
    $hash{$type}{$pre}++;
}else{
    my @ids = split(/,/,$id);
    foreach my $id2 (@ids){
      $type = ${${$all_genes{$id2}}[0]}{"gene_type"} if($id2 ne '*' and $id2 ne 'NA');
      die "$_\n" if(not defined $type);
      $hash{$type}{$pre}++;
    }
}
    $hash{$assign}{$pre}++;
    $hash{MultiOverlapping}{$pre}++ if($num > 1 and $assign ne 'Unassigned_Secondary');
    $seen{$read} = 1;
    $seen{$read} = 2 if($id ne '*' and $id ne 'NA' and $assign ne 'Unassigned_Secondary');
  }
  $hash{TotalReads}{$pre}=scalar keys %seen;
  $hash{ReadsAssigned}{$pre} = scalar( grep { $seen{$_} == 2 } keys %seen );
  $hash{ReadsUnassigned}{$pre} = scalar( grep { $seen{$_} == 1 } keys %seen );
  close $in;
}
closedir($din);
print "Assignment";
foreach my $pre (sort @samp){
  print "\t$pre";
}
print "\n";
foreach my $t (sort keys %hash){
  print "$t";
  foreach my $pre (sort @samp){
    if(defined $hash{$t}{$pre}){
      print "\t$hash{$t}{$pre}";
    }else{
      print "\t0";
    }
  }
  print "\n";
}
exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
# USAGE:
# perl make_rpkm_table.pl -t /nfs/bio/db/Homo_sapien/gencode19_with_contaminates/avgTranscriptLength.tsv -m /nfs/labs/nephew/mcf7_timecourse/mcf7_valid.final.tsv -c /nfs/labs/nephew/mcf7_timecourse/summary/mcf7_timecourse_Hist.antisense.gene.cnts -o mcf7_timecourse.rpkm.tsv

my %options;
my $results = GetOptions(\%options, 
                          "transLength|t=s", # Avg Transcript Length
                          "counts|c=s",      # Counts Files Per Transcript
                          "output|o=s",      # Output File
                          "mappedCount|m=s", # Total Mapped (Valid) Reads
                        );
open(my $tin, $options{transLength}) or die "can't open $options{transLength}: $!\n";
my %transLength;
while(<$tin>){
#  next unless($_ =~ /^ENSG/o);
  chomp $_;
  my @line = split(/\t/, $_);
  $transLength{$line[0]} = $line[1];
}
close $tin;

open(my $min, "$options{mappedCount}") or die "Could not open $options{mappedCount}: $!\n";
my %mappedCount;
while(<$min>){
  chomp $_;
  my @line = split(/\t/, $_);
  $mappedCount{$line[0]} = $line[1];
}
close($min);

open(my $cin, "$options{counts}") or die "can't open $options{counts}: $!\n";
open(my $out, ">$options{output}") or die "can't open $options{output}: $!\n";
my $header = <$cin>;
chomp $header;
my @h = split(/\t/, $header);
my %headerLocMap;
for(my $i=0; $i <= scalar @h; $i++){
  $headerLocMap{$i} = $h[$i];
}
print $out $header, "\n";
while (<$cin>){
  chomp $_;
  next unless($_ =~ /^ENS/);
#  next unless($_ =~ /^ENSMUSG/);
  my @line = split /\t/, $_;
  my $id = $line[0];
  print $out $id;
  my $translength = $transLength{$id};
  for(my $i=1; $i < scalar @line; $i++){
    my $genecount = $line[$i];

    ### RPKM calculation for a gene in each sample iteratively
    # RPKM=(number of mapping reads)*1000bp*1 million reads/[(length of transcript)*(number of total reads )]
    my $totalreads = $mappedCount{$headerLocMap{$i}};
    my $num = 1000000000 * $genecount;
    my $den = $translength * $totalreads;
#    my $rpkm = (($num / $den) + 1);  ### psuedocounts added
    my $rpkm = $num / $den;  ### psuedocounts not added
    print $out "\t$rpkm";
  }
  print $out "\n";
}
close $cin;
close $out;
exit;

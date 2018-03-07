#!/usr/bin/env perl
use strict;
use warnings;

my $dir = shift;
#my $map = shift;
#open(my $min, $map) or die "can't open $map: $!\n";
#my %labels;
#while(<$min>){
#  chomp $_;
#  my @line = split(/\t/, $_);
#  $labels{$line[0]} = $line[1];
#}
#close $min;

opendir(my $din, $dir) or die "$!\n";
while(my $file = readdir($din)){
## {
##   "label" : "bam_alignments",
##   "key" : "BAM alignments",
##   "storeClass" : "JBrowse/Store/SeqFeature/BAM",
##   "urlTemplate" : "../../simulated-sorted.bam",
##   "type" : "Alignments2"
## }
## MCF7_Z_vector_smRNA.3.sorted.bam
  next unless($file =~ /\.sorted\.bam$/o);
  $file =~ /(\S+)\.sorted\.bam/o;
  my $label = $1;
  my $bam = $dir . '/'. $file;
  my $bai = $dir . '/'. $file . '.bai';

  print '      {', "\n";
  print '        "label" : "'.$label.' BAM alignments",', "\n";
  print '        "key" : "'.$label.' BAM alignments",', "\n";
  print '        "category" : "BAM Files",',"\n";
  print '        "storeClass" : "JBrowse/Store/SeqFeature/BAM",',"\n";
  print '        "urlTemplate" : "bam/'.$file, "\",\n";
#  print '        "baiUrlTemplate" : "'.$bai, "\",\n\n";
  print '        "type" : "Alignments2"', "\n";
  print '      },',"\n";

}
close $din;
exit;

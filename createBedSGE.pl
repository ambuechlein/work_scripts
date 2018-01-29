#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "directory|d=s",
                         "outdirectory|o=s",
                        );


opendir(my $din, $options{directory}) or die "can't open directory: $!\n";
while(my $sub=readdir($din)){
  next unless(-d "$dir/$sub" and $sub =~ /_thout$/o);
  my $lib = $sub;
  $lib =~ s/\.fastq_thout//go;

  my $bam =  "$options{directory}/$sub/accepted_hits_sorted.bam";
  my $bed =  "$options{outdirectory}/$lib.bed";
  my $nbed = "$options{outdirectory}/$lib.negative.bed";
  my $pbed = "$options{outdirectory}/$lib.positive.bed";

  open(my $sh, ">sam4bed.$lib.sh") or die "can't open sge script for $lib: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo "Creating BED"', "\n";
  print $sh "/home/abuechle/bin/bedtools-2.17.0/bin/bamToBed -split -i $bam > $bed" .' || { echo "Bam2Bed failed"; exit 1; }'."\n";
  print $sh 'echo "make negative BED"', "\n";
  print $sh "perl -nle 'chomp \$_; print \$_ if(\$_ =~ /\\-\$/o);' $bed > $nbed" .' || { echo "negative bed creation failed"; exit 1; }'."\n";
  print $sh 'echo "make positive BED"', "\n";
  print $sh "perl -nle 'chomp \$_; print \$_ if(\$_ =~ /\\+\$/o);' $bed > $pbed" .' || { echo "positive bed creation failed"; exit 1; }'."\n";
  print $sh "rm $bed\n";
  print $sh 'echo "make negative BEDGraph"', "\n";
  print $sh "/home/abuechle/bin/bedtools-2.17.0/bin/genomeCoverageBed -bg -split -trackline -trackopts 'name=\"$lib negative\" visibility=2 color=255,30,30' -i $nbed -g /home/abuechle/Desktop/hg19.chrom.sizes > $nbed"."graph" .' || { echo "negative bedgraph creation failed"; exit 1; }'."\n";
  print $sh 'echo "make positive BEDGraph"', "\n";
  print $sh "/home/abuechle/bin/bedtools-2.17.0/bin/genomeCoverageBed -bg -split -trackline -trackopts 'name=\"$lib positive\" visibility=2 color=255,30,30' -i $pbed -g /home/abuechle/Desktop/hg19.chrom.sizes > $pbed"."graph" .' || { echo "positive bedgraph creation failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 sam4bed.$lib.sh");
#  system("qsub -q cluster@ccf17,cluster@ccf18 sam4bed.$lib.sh");
}
close $in;
exit;

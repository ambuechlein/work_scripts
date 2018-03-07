#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "outdir|o=s",
                        "sampleorder|s=s",
#                        "chrSizes|c=s",
                       );

my %map;
open(my $sin, $options{sampleorder}) or die "Can't open $options{sampleorder}: $!\n";
while (<$sin>) {
  chomp $_;
  my @d = split(/\t/, $_);
  $map{$d[0]} = $d[1];
}
close $sin;

opendir(my $din, $options{dir}) or die "can't open input directory $options{dir}: $!\n";
system("mkdir $options{outdir}") unless(-e $options{outdir});
while(my $dir = readdir($din)){
  next unless( $dir =~ /\.paired\.fastq\.gz_hisat$/o);
  $dir =~ /(\S+)\.paired\.fastq\.gz_hisat/o;
  die "couldn't parse basename correctly for $dir\t$1\n" unless(defined $map{$1});
  my $pre = $1;
  my $base = $map{$1};
  $base =~ s/(\(|\))/_/go;
  $base =~ s/_+/_/go;

  my $bam = "$options{dir}/$dir/accepted_hits.sorted.bam";
  my $nbw =   $options{outdir} . "/$base.negative.bw";
  my $pbw =   $options{outdir} . "/$base.positive.bw";

  open(my $sh, ">bam2bw.$pre.sh") or die "can't open sge script for $base: $!\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> Bam2BW "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "make positive bigwig"', "\n";
  print $sh "bamCoverage -b $bam --filterRNAstrand forward -p 8 -of bigwig -o $pbw" .' || { echo "positive bigwig creation failed"; exit 1; }'."\n";
  print $sh 'echo "make negative bigwig"', "\n";
  print $sh "bamCoverage -b $bam --filterRNAstrand reverse -p 8 -of bigwig -o $nbw" .' || { echo "negative bigwig creation failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 bam2bw.$pre.sh");
#  system("qsub -q bigmem -pe pe_slots 8 -wd $cwd bam2bw.$pre.sh");
}

closedir $din;
exit;

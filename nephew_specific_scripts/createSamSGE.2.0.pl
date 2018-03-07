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
                        "chrSizes|c=s",
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
while(my $sam = readdir($din)){
  next unless( $sam =~ /\.sam.gz$/o);
  $sam =~ /(\S+)\.sam.gz/o;
  die "couldn't parse basename correctly for $sam\t$1\n" unless(defined $map{$1});
  my $pre = $1;
  my $base = $map{$1};
  $base =~ s/(\(|\))/_/go;
  $base =~ s/_+/_/go;

  $sam = $options{dir} . '/' . $sam;
  my $samout = $options{outdir} . '/' . $base . '.fixed.sam';
  my $bamout = $options{outdir} . '/' . $base . '.bam';
  my $sort =   $options{outdir} . '/' . $base . '.sorted';
  my $bed =    $options{outdir} . "/$base.bed";
  my $nbed =   $options{outdir} . "/$base.negative.bed";
  my $pbed =   $options{outdir} . "/$base.positive.bed";
  
  my $command1 = "samtools view -bS  -@ 8 $sam > $bamout";
#  my $command2 = "samtools sort -@ 8 -m 1G $bamout $sort";
  my $command2 = "samtools sort -@ 8 -m 1G -O BAM -o $sort.bam $bamout";
  my $command3 = "samtools index $sort.bam";

  open(my $sh, ">sam.$pre.sh") or die "can't open sge script for $base: $!\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> Sam4Bedgraph "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";


  print $sh 'echo "Starting Sam 2 bam"', "\n";
  print $sh $command1 . ' || { echo "samtools 2 bam failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting Bam sort"', "\n";
  print $sh $command2 . ' || { echo "samtools sort failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting Bam index"', "\n";
  print $sh $command3 . ' || { echo "samtools index failed"; exit 1; }' . "\n";

  print $sh 'echo "make negative BEDGraph"', "\n";
  print $sh "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -split -strand - -trackline -trackopts 'name=\"$base negative\" visibility=2 color=255,30,30' -ibam $sort.bam -g $options{chrSizes} > $nbed"."graph" .' || { echo "negative bedgraph creation failed"; exit 1; }'."\n";
  print $sh 'echo "make positive BEDGraph"', "\n";
  print $sh "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -split -strand + -trackline -trackopts 'name=\"$base positive\" visibility=2 color=255,30,30' -ibam $sort.bam -g $options{chrSizes} > $pbed"."graph" .' || { echo "positive bedgraph creation failed"; exit 1; }'."\n";
  print $sh "rm $bamout\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 sam.$pre.sh");
  system("qsub -q cluster -pe pe_slots 8 -wd $cwd sam.$pre.sh");
}

closedir $din;
exit;

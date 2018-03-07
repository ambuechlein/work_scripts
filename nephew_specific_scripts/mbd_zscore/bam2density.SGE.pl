#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "outdir|o=s",
                       );

opendir(my $din, $options{dir}) or die "can't open input directory $options{dir}: $!\n";
while(my $bam = readdir($din)){
  next unless( $bam =~ /\.bam$/o);
  $bam =~ /(\S+)\.bam/o;
  my $pre = $1;

  $bam = $options{dir} . '/' . $bam;
  my $bed = $options{outdir} . "/$pre.bed";
  my $den = $options{outdir} . "/$pre.density";
  
  my $command = "bamToBed -split -i $bam | perl /nfs/labs/nephew/scripts/mbd_zscore/newDensityByChr.pl -i - -o $den -s 1";

  open(my $sh, ">bam2density.$pre.sh") or die "can't open sge script for $pre: $!\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> Bam2Density "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh "mkdir -p $options{outdir}/byChr\n";
  print $sh 'echo "starting bam2density"', "\n";
  print $sh $command . ' || { echo "bam2density failed"; exit 1; }' . "\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 bam2density.$pre.sh");
  system("qsub -q bigmem\@intron,bigmem\@antigen -pe pe_slots 8 -wd $cwd bam2density.$pre.sh");
}

closedir $din;
exit;

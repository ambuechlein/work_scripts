#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "out|o=s",
                       );
my $dir = $options{dir};
my $outdir = $options{out};
my $cwd = cwd();

opendir(my $din, $dir) or die "can't open directroy: $!\n";
while(my $f = readdir($din)){
  next unless($f =~ /.fastq.gz$/o);
  my $prefix = $f;  $prefix =~ s/\.fastq.gz$//o;
  my $file = $dir . '/' . $f;
  my $fullout = $outdir.'/'.$prefix.'/';

  my $command9 = "/home/abuechle/bin/bedtools2/bin/coverageBed -counts -F .2 -a /nfs/labs/nephew/human_databases/gencode_v23/500bpwindows.bed -b ${fullout}${prefix}.chr.sorted.bam > ${fullout}${prefix}.windowCounts.tsv";

  open(my $sh, ">countwindow.$f.sh") or die "can't open shell script countwindow.$f.sh: $!\n\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> count "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting window count"', "\n";
  print $sh $command9 . ' || { echo "count failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 countwindow.$f.sh");
  system("qsub -q bigmem\@intron,bigmem\@exon.cgb.indiana.edu -pe pe_slots 8 -wd $cwd countwindow.$f.sh");
}
exit;

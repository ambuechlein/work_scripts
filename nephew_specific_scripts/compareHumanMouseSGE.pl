#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my %options;
my $results = GetOptions(\%options,
                         "human|h=s",
                         "mouse|m=s",
                         "out|o=s",
                        );

my $cwd = cwd();
#*.chr.sorted.bam
opendir(my $hd, $options{human}) or die "Can't open directory $options{human}: $!\n";
while(my $pre = readdir($hd)){
  next unless($pre =~ /^GSF/o); 
  next unless(-d $pre);
  my $hbam = "$options{human}/$pre/$pre.sorted.bam";
  die "$hbam does not exist\n" unless(-e $hbam);
  my $mbam = "$options{mouse}/$pre/$pre.sorted.bam";
  die "$mbam does not exist\n" unless(-e $mbam);
  my $command = "perl /nfs/labs/nephew/Pili/bowtie/compareHumanMouse.pl -h $hbam -m $mbam -o $options{out}";

  open(my $sh, ">humanMouseCounts.$pre.sh") or die "can't open shell script humanMouseCounts.$pre.sh: $!\n\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> countHumanMouse "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh "$command" . ' || { echo "counts failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 humanMouseCounts.$pre.sh");
  system("qsub -q bigmem\@intron,bigmem\@exon.cgb.indiana.edu -pe pe_slots 4 -wd $cwd humanMouseCounts.$pre.sh");
}
closedir $hd;

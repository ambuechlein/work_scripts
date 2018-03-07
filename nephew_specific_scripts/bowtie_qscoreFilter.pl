#!/usr/bin/env perl
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

  my $command1 = "bowtie2 -p 8 -k 20 -S ${fullout}${prefix}.sam --un ${fullout}${prefix}.notmapped -x /nfs/labs/nephew/human_databases/gencode_v23/genome -U $file";
  my $command2 = "samtools view -q 40 -@ 8 -bS -o ${fullout}${prefix}.bam ${fullout}${prefix}.sam";
  my $command3 = "samtools view -h -F4 ${fullout}${prefix}.bam | awk '{if(match(\$0, /^\\\@/) || match(\$3, /^chr/) ){print \$0}}' | samtools view -Sb - > ${fullout}${prefix}.chr.bam";
  my $command4 = "samtools sort -@ 8 -m 1G -O bam -o ${fullout}${prefix}.sorted.bam ${fullout}${prefix}.chr.bam";
  my $command5 = "samtools index ${fullout}${prefix}.sorted.bam";
  my $command6 = "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -trackline -trackopts 'name=\"$prefix\" visibility=2 color=255,30,30' -ibam ${fullout}${prefix}.sorted.bam -g /nfs/labs/nephew/human_databases/gencode_v23/chrSizes.tsv > ${fullout}${prefix}.bedgraph";
  my $command7 = "/home/abuechle/bin/bedtools2/bin/coverageBed -counts -F .2 -a /nfs/labs/nephew/human_databases/gencode_v25/500bpwindows.bed -b ${fullout}${prefix}.sorted.bam > ${fullout}${prefix}.windowCounts.tsv";


  open(my $sh, ">bowtie.$f.sh") or die "can't open shell script bowtie.$f.sh: $!\n\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> bowtie2 "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'mkdir -p '.$fullout, "\n";
  print $sh 'echo "Starting Bowtie"', "\n";
  print $sh $command1 . ' || { echo "bowtie failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Qscore Filter"', "\n";
  print $sh $command2 . ' || { echo "best filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting chr Filter"', "\n";
  print $sh $command3 . ' || { echo "samtools chr filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools sort"', "\n";
  print $sh $command4 . ' || { echo "samtools sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools index 1"', "\n";
  print $sh $command5 . ' || { echo "samtools index failed"; exit 1; }'."\n";
  print $sh 'echo "Starting bedgraph"', "\n";
  print $sh $command6 . ' || { echo "bedgraph failed"; exit 1; }'."\n";
  print $sh 'echo "Starting window count"', "\n";
  print $sh $command7 . ' || { echo "window count failed"; exit 1; }'."\n";
  print $sh 'echo "Starting compressing and deleting"', "\n";
  print $sh "pigz -f ${fullout}${prefix}.sam ${fullout}${prefix}.notmapped" . ' || { echo "gzip failed"; exit 1; }'."\n";
  print $sh "rm -f ${fullout}${prefix}.bam ${fullout}${prefix}.chr.bam" . ' || { echo "rm failed"; exit 1; }'."\n";

  print $sh 'echo "Starting summary counts"', "\n";
  print $sh "/nfs/labs/nephew/countMappedReads.pl ${fullout}${prefix}.sorted.bam" . ' || { echo "summary count failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 bowtie.$f.sh");
  system("qsub -q bigmem -pe pe_slots 4 -wd $cwd bowtie.$f.sh");
}
exit;

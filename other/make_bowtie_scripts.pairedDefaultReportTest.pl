#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "genome|g=s",
                        "out|o=s",
                       );
my $dir = $options{dir};
my $outdir = $options{out};

opendir(my $din, $dir) or die "can't open directroy: $!\n";
while(my $f = readdir($din)){
  next unless($f =~ /R1_001.paired.fastq.gz$/o);
  my $file1 = "$dir/$f";
  my $file2 = $file1;
  $file2 =~ s/R1_001.paired.fastq.gz/R2_001.paired.fastq.gz/go;
  my $prefix = $f;
  $prefix =~ s/_R1_001.paired.fastq.gz//g;
  my $fullout = $outdir.'/'.$prefix.'/';

  my $command1 = "/N/u/abuechle/Karst/miniconda3/bin/bowtie2 -p 8 --mm --very-sensitive -S ${fullout}${prefix}.default.sam -x $options{genome} -1 $file1 -2 $file2";
  my $command2 = "samtools view -@ 8 -bS -F 256 -o ${fullout}${prefix}.default.bam ${fullout}${prefix}.default.sam";
#  my $command3 = "perl /N/u/abuechle/Karst/bin/filterSam.pl ${fullout}${prefix}.bam | samtools view -h -bS -o ${fullout}${prefix}.best.bam - ";
  my $command4 = "samtools sort -@ 8 -m 4G -O bam -o ${fullout}${prefix}.sorted.bam ${fullout}${prefix}.default.bam";
  my $command5 = "samtools index ${fullout}${prefix}.sorted.bam";
#   my $command6 = " /N/u/abuechle/Karst/bin/bedtools2-2.26.0/bin/genomeCoverageBed -bg -ibam ${fullout}${prefix}.sorted.bam -g /N/dc2/projects/cgbgsf/Rnor_6.0/chrSizes.tsv > ${fullout}${prefix}.bedgraph";
#   my $command7 = "perl /N/u/abuechle/Karst/bin/singleMappedBam.random.pl ${fullout}${prefix}.sorted.bam | samtools view -@ 8 -bS - > ${fullout}${prefix}.s.bam";
#   my $command8 = "samtools sort -@ 8 -m 4G -O bam -o ${fullout}${prefix}.s.sorted.bam ${fullout}${prefix}.s.bam";
#   my $command9 = "/N/u/abuechle/Karst/bin/bedtools2-2.26.0/bin/genomeCoverageBed -bg -ibam ${fullout}${prefix}.s.sorted.bam -g /N/dc2/projects/cgbgsf/Rnor_6.0/chrSizes.tsv > ${fullout}${prefix}.single.bedgraph";

  open(my $sh, ">bowtieDF.$f.sh") or die "can't open shell script bowtieDF.$f.sh: $!\n\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> bowtie2 "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'mkdir -p '.$fullout, "\n";
  print $sh 'echo "Starting Bowtie"', "\n";
  print $sh $command1 . ' || { echo "bowtie failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Primary Filter"', "\n";
  print $sh $command2 . ' || { echo "Primary Map Filiter failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting Best Filter"', "\n";
#  print $sh $command3 . ' || { echo "best filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools sort 1"', "\n";
  print $sh $command4 . ' || { echo "samtools sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools index 1"', "\n";
  print $sh $command5 . ' || { echo "samtools index failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting bam2bedgraph 1"', "\n";
#  print $sh $command6 . ' || { echo "bam2bedgraph1 failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting single mapped"', "\n";
#  print $sh $command7 . ' || { echo "single mapped failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting single map sort"', "\n";
#  print $sh $command8 . ' || { echo "single sort failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting single bam2bed"', "\n";
#  print $sh $command9 . ' || { echo "single bam2bed failed"; exit 1; }'."\n";

#  print $sh 'echo "Starting bedgraph sort"', "\n";
#  print $sh "sort --buffer-size=50G -k1,1 -k2,2n ${fullout}${prefix}.bedgraph > ${fullout}${prefix}.bedgraph2" . ' || { echo "bedgraph sort1 failed"; exit 1; }'."\n";
#  print $sh "mv ${fullout}${prefix}.bedgraph2 ${fullout}${prefix}.bedgraph" . ' || { echo "bedgraph move1 failed"; exit 1; }'."\n";
#  print $sh "sort --buffer-size=50G -k1,1 -k2,2n ${fullout}${prefix}.single.bedgraph > ${fullout}${prefix}.single.bedgraph2" . ' || { echo "bedgraph sort1 failed"; exit 1; }'."\n";
#  print $sh "mv ${fullout}${prefix}.single.bedgraph2 ${fullout}${prefix}.single.bedgraph" . ' || { echo "bedgraph move1 failed"; exit 1; }'."\n";

  print $sh 'echo "Starting compressing and deleting"', "\n";
#  print $sh "pigz -f ${fullout}${prefix}.sam" . ' || { echo "gzip failed"; exit 1; }'."\n";
  print $sh "rm -f ${fullout}${prefix}.default.sam ${fullout}${prefix}.default.bam" . ' || { echo "rm failed"; exit 1; }'."\n";

  print $sh 'echo "Starting summary counts"', "\n";
  print $sh "/N/u/abuechle/Karst/bin/work_scripts/other/countMappedReads.pl ${fullout}${prefix}.sorted.bam" . ' || { echo "summary count failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 bowtieDF.$f.sh");
}
exit;

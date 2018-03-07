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

  my $command1 = "bowtie2 -p 8 -k 20 -S ${fullout}${prefix}.sam --un ${fullout}${prefix}.notmapped -x /nfs/labs/nephew/human_databases/gencode_m12/genome -U $file";
  my $command2 = "perl /nfs/labs/nephew/filterSam.pl ${fullout}${prefix}.sam > ${fullout}${prefix}.best.sam";

  my $command3 = "samtools view -@ 8 -bS ${fullout}${prefix}.best.sam > ${fullout}${prefix}.bam";
  my $command4 = "samtools sort -@ 8 -m 1G -O bam -o ${fullout}${prefix}.sorted.bam ${fullout}${prefix}.bam";
  my $command5 = "samtools index ${fullout}${prefix}.sorted.bam";
  my $command6 = "samtools view -h -F4 ${fullout}${prefix}.sorted.bam | awk '{if(match(\$0, /^\\\@/) || match(\$3, /^chr/) ){print \$0}}' | samtools view -Sb - > ${fullout}${prefix}.chr.bam";
  my $command7 = "samtools sort -@ 8 -m 1G -O bam -o ${fullout}${prefix}.chr.sorted.bam ${fullout}${prefix}.chr.bam";

#  my $command8 = "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -trackline -trackopts 'name=\"$prefix\" visibility=2 color=255,30,30' -ibam ${fullout}${prefix}.chr.sorted.bam -g /nfs/labs/nephew/human_databases/gencode_m12/chrSizes.tsv > ${fullout}${prefix}.chr.bedgraph";
  my $command8 = "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -ibam ${fullout}${prefix}.chr.sorted.bam -g /nfs/labs/nephew/human_databases/gencode_m12/chrSizes.tsv > ${fullout}${prefix}.chr.bedgraph";


  my $command10 = "perl /nfs/labs/nephew/scripts/singleMappedBam.pl ${fullout}${prefix}.sorted.bam | samtools view -@ 8 -bS - > ${fullout}${prefix}.s.bam";
  my $command11 = "samtools sort -@ 8 -m 1G -O bam -o ${fullout}${prefix}.s.sorted.bam ${fullout}${prefix}.s.bam";
# my $command12 = "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -trackline -trackopts 'name=\"$prefix\" visibility=2 color=255,30,30' -ibam ${fullout}${prefix}.s.sorted.bam -g /nfs/labs/nephew/human_databases/gencode_m12/chrSizes.tsv > ${fullout}${prefix}.single.bedgraph";
  my $command12 = "/home/abuechle/bin/bedtools2/bin/genomeCoverageBed -bg -ibam ${fullout}${prefix}.s.sorted.bam -g /nfs/labs/nephew/human_databases/gencode_m12/chrSizes.tsv > ${fullout}${prefix}.single.bedgraph";

  my $command13 = "/home/abuechle/bin/bedtools2/bin/coverageBed -counts -F .2 -a /nfs/labs/nephew/human_databases/gencode_m12/500bpwindows.bed -b ${fullout}${prefix}.s.sorted.bam > ${fullout}${prefix}.s.windowCounts.tsv";
  my $command14 = "/home/abuechle/bin/bedtools2/bin/coverageBed -counts -F .2 -a /nfs/labs/nephew/human_databases/gencode_m12/500bpwindows.bed -b ${fullout}${prefix}.chr.sorted.bam > ${fullout}${prefix}.windowCounts.tsv";


  open(my $sh, ">bowtie.$f.sh") or die "can't open shell script bowtie.$f.sh: $!\n\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> bowtie2 "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh "which bowtie2\n";
  print $sh 'mkdir -p '.$fullout, "\n";
  print $sh 'echo "Starting Bowtie"', "\n";
  print $sh $command1 . ' || { echo "bowtie failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Filter"', "\n";
  print $sh $command2 . ' || { echo "best filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools conversions"', "\n";
  print $sh $command3 . ' || { echo "samtools sam2bam failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools sort 1"', "\n";
  print $sh $command4 . ' || { echo "samtools sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools index 1"', "\n";
  print $sh $command5 . ' || { echo "samtools index failed"; exit 1; }'."\n";
  print $sh 'echo "Starting awk"', "\n";
  print $sh $command6 . ' || { echo "awk filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools sort 2"', "\n";
  print $sh $command7 . ' || { echo "sort filter failed"; exit 1; }'."\n";
  print $sh 'echo "Starting samtools index 2"', "\n";
  print $sh "samtools index ${fullout}${prefix}.chr.sorted.bam" . ' || { echo "samtools index2 failed"; exit 1; }'."\n";
  print $sh 'echo "Starting bedgraph conversions"', "\n";
  print $sh $command8 . ' || { echo "bed2bedgraph failed"; exit 1; }'."\n";
  print $sh 'echo "Starting single mapped"', "\n";
  print $sh $command10 . ' || { echo "single mapped failed"; exit 1; }'."\n";
  print $sh 'echo "Starting single sort"', "\n";
  print $sh $command11 . ' || { echo "single sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting single bedgraph"', "\n";
  print $sh $command12 . ' || { echo "single window count failed"; exit 1; }'."\n";
  print $sh 'echo "Starting window Count"', "\n";
  print $sh $command13 . ' || { echo "single mapped bedgraph failed"; exit 1; }'."\n";
  print $sh 'echo "Starting window Count2"', "\n";
  print $sh $command14 . ' || { echo "window count failed"; exit 1; }'."\n";

  print $sh 'echo "Starting bedgraph sort"', "\n";
  print $sh "sort -k1,1 -k2,2n ${fullout}${prefix}.chr.bedgraph > ${fullout}${prefix}.chr.bedgraph2" . ' || { echo "bedgraph sort1 failed"; exit 1; }'."\n";
  print $sh "mv ${fullout}${prefix}.chr.bedgraph2 ${fullout}${prefix}.chr.bedgraph" . ' || { echo "bedgraph move1 failed"; exit 1; }'."\n";
  print $sh "sort -k1,1 -k2,2n ${fullout}${prefix}.single.bedgraph > ${fullout}${prefix}.single.bedgraph2" . ' || { echo "bedgraph sort1 failed"; exit 1; }'."\n";
  print $sh "mv ${fullout}${prefix}.single.bedgraph2 ${fullout}${prefix}.single.bedgraph" . ' || { echo "bedgraph move1 failed"; exit 1; }'."\n";


  print $sh 'echo "Starting compressing and deleting"', "\n";
  print $sh "pigz -f ${fullout}${prefix}.sam" . ' || { echo "gzip failed"; exit 1; }'."\n";
  print $sh "rm -f ${fullout}${prefix}.bam ${fullout}${prefix}.best.sam ${fullout}${prefix}.chr.bam ${fullout}${prefix}.s.bam" . ' || { echo "rm failed"; exit 1; }'."\n";

  print $sh 'echo "Starting summary counts"', "\n";
  print $sh "/nfs/labs/nephew/countMappedReads.pl ${fullout}${prefix}.sorted.bam" . ' || { echo "summary count failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 bowtie.$f.sh");
  system("qsub -q bigmem -pe pe_slots 4 -wd $cwd bowtie.$f.sh");
}
exit;

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
while(my $in = readdir($din)){
# 130604_Pool3_1_ATCACG_L003_R1.fixed.sam.gz
  next unless( $in =~ /.sorted.bam$/o);
  $in =~ /(\S+)\.sorted.bam/o;
  my $pre = $1;

  my $bamin = "$options{dir}/$in";
  my $samout = $options{outdir} . '/' . $pre . '.sam';
  my $bamout = $options{outdir} . '/' . $pre . '.bam';
  my $sort =   $options{outdir} . '/' . $pre . '.sorted';
  
  my $command1 = "/home/abuechle/bin/sizeFilterBam.pl -i $bamin -o $samout";
  my $command6 = "samtools view -@ 16 -bS -h $samout > $bamout";
  my $command2 = "samtools sort -@ 16 -m 4G $bamout $sort";
  my $command3 = "samtools index $sort.bam";
  my $command4 = "/home/abuechle/bin/samtools-0.1.19/samtools rmdup -s $sort.bam $options{outdir}/$pre.deduped.bam";
  my $command5 = "samtools sort -@ 16 -m 4G $options{outdir}/$pre.deduped.bam $options{outdir}/$pre.deduped.sorted";
  my $command7 = "htseq-count -f bam $sort.bam /nfs/bio/db/Homo_sapien/gencode_v19/gencode.v19.annotation.gtf > $options{outdir}/$pre.htseqcount.tsv";
  my $command8 = "htseq-count -f bam $options{outdir}/$pre.deduped.sorted.bam /nfs/bio/db/Homo_sapien/gencode_v19/gencode.v19.annotation.gtf > $options{outdir}/$pre.deduped.htseqcount.tsv";

  open(my $sh, ">sam.$pre.sh") or die "can't open sge script for $pre: $!\n";
  print $sh '#!/bin/bash', "\n";
  print $sh 'echo ">>>>> Sam4Bedgraph "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting size filter 2 sam"', "\n";
  print $sh $command1 . ' || { echo "bam size filter failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting sam 2 bam"', "\n";
  print $sh $command6 . ' || { echo "samtools sam 2 bam failed" exit 1; }' . "\n";
  print $sh 'echo "Starting Bam sort"', "\n";
  print $sh $command2 . ' || { echo "samtools sort failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting Bam index"', "\n";
  print $sh $command3 . ' || { echo "samtools index failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting rmdup"', "\n";
  print $sh $command4 . ' || { echo "samtools rmdup failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting dedupe sort"', "\n";
  print $sh $command5 . ' || { echo "samtools sort rmdup failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting deduped index"', "\n";
  print $sh "samtools index $options{outdir}/$pre.deduped.sorted.bam" . ' || { echo "samtools deduped index failed"; exit 1; }' . "\n";
  print $sh 'echo "Removing tmp files"', "\n";
  print $sh "rm $samout $bamout $options{outdir}/$pre.deduped.bam\n";
  print $sh 'echo "Starting htseq"', "\n";
  print $sh $command7 . ' || { echo "htseq failed"; exit 1; }' . "\n";
  print $sh 'echo "Starting htseq dedupe"', "\n";
  print $sh $command8 . ' || { echo "htseq dedupe failed"; exit 1; }' . "\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 sam.$pre.sh");
  system("qsub -q bigmem\@intron,bigmem\@exon.cgb.indiana.edu -pe pe_slots 8 -wd $cwd sam.$pre.sh");
}

closedir $din;
exit;

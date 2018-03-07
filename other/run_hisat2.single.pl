#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "out|o=s",
                        "genome|g=s",
                       );
my $dir = $options{dir};
my $out = $options{out};

opendir(my $din, $dir) or die "can't open directroy: $!\n";
while(my $f = readdir($din)){
  next unless($f =~ /fastq.gz$/o);
  my $file1 = "$dir/$f";
  my $outdir = $out.'/'.$f.'_hisat';
  my $pre = $f;
  $pre =~ s/_R1_001.fastq.gz//g;
  my $command = "hisat2 -q --rna-strandness R --new-summary --summary-file $outdir/$pre.summary.tsv --un-gz $outdir/unmapped.tsv.gz -p 8 --mm -x $options{genome} -U $file1 -S $outdir/accepted_hits.sam";

  open(my $sh, ">hisat.$f.sh") or die "can't open shell script $f: $!\n\n";

  print $sh '#!/bin/bash';
  print $sh "\n";
  print $sh 'echo ">>>>>  Hisat"', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'mkdir -p '.$outdir . ' || { echo "mkdir failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Hisat"', "\n";
  print $sh $command . ' || { echo "hisat failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Sam2Bam"', "\n";
  print $sh "samtools view -bS -@ 8 $outdir/accepted_hits.sam > $outdir/accepted_hits.bam" . ' || { echo "sam2bamfailed"; exit 1; }'."\n";
  print $sh 'echo "Starting Sam Sort"', "\n";
  print $sh "samtools sort -@ 8 -m 1G -O BAM -o $outdir/accepted_hits.sorted.bam $outdir/accepted_hits.bam" . ' || { echo "sam sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Sam Index"', "\n";
  print $sh "samtools index $outdir/accepted_hits.sorted.bam" . ' || { echo "sam sort failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Delete"', "\n";
  print $sh "rm $outdir/accepted_hits.sam $outdir/accepted_hits.bam" . ' || { echo "delete failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 hisat.$f.sh");
#  system("qsub -q bigmem -pe pe_slots 8 -wd $cwd hisat.$f.sh");
}
exit;

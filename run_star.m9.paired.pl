#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "out|o=s",
                       );
my $dir = $options{dir};
my $out = $options{out};

opendir(my $din, $dir) or die "can't open directroy: $!\n";
while(my $f = readdir($din)){
# GSF1135-1_S1_R1_001.paired.fastq.gz
# GSF1135-1_S1_R1_001.unpaired.fastq.gz
# GSF1135-1_S1_R2_001.paired.fastq.gz
# GSF1135-1_S1_R2_001.unpaired.fastq.gz
  next unless($f =~ /R1_001.paired.fastq.gz/o);
  next unless($f =~ /\.fastq.gz$/o);
  my $f2 = $f;
#  $f2 =~ s/_R1//go;
  $f2 =~ s/paired/unpaired/go;
  my $f3 = $f;
  $f3 =~ s/_R1_/_R2_/go;
  my $f4 = $f3;
  $f4 =~ s/paired/unpaired/go;
  my $file = "$dir/$f $dir/$f3";
  my $outdir = $out.'/'.$f.'_star';
  my $pre = $f;
  $pre =~ s/_001\S+\.fastq.gz$/_001/;

  my $command = "/home/abuechle/bin/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir /nfs/labs/nephew/human_databases/gencode_m9/ --readFilesIn $file --readFilesCommand zcat --outFilterMultimapNmax 25 --outFileNamePrefix $outdir/ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts";
  open(my $sh, ">star.$f.sh") or die "can't open shell script $f: $!\n\n";

  print $sh '#!/bin/bash';
  print $sh "\n";
  print $sh 'echo ">>>>> Star "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'mkdir -p '.$outdir . ' || { echo "mkdir failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Star"', "\n";
  print $sh $command . ' || { echo "star failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting Samtools"', "\n";
#  print $sh '/nfs/bio/sw/bin/samtools sort -@ 8 -m 1G ' .$outdir. '/accepted_hits.bam ' . $outdir . '/accepted_hits.sorted' . ' || { echo "samtools sort failed"; exit 1; }'."\n";
#  print $sh '/nfs/bio/sw/bin/samtools index '.$outdir. '/accepted_hits.sorted.bam' . ' || { echo "samtools index failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 star.$f.sh");
  system("qsub -q bigmem,cluster -pe pe_slots 8 -wd $cwd star.$f.sh");
}
exit;

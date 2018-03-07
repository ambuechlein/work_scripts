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
# GSF1016-polyA-2717_S7_R1_001.fastq.gz
# GSF1016-polyA-2746_S1_R1_001.fastq.gz
# GSF1016-polyA-2747_S8_R1_001.fastq.gz
# GSF1016-polyA-3591_S4_R1_001.fastq.gz
  next if($f =~ /small/o);
  next unless($f =~ /ishi/o);
  next unless($f =~ /\.fastq.gz$/o);
  my $file = "$dir/$f";
  my $outdir = $out.'/'.$f.'_star';
  my $pre = $f;
  $pre =~ s/\.fastq.gz$//;

  my $command = "/home/abuechle/bin/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --outReadsUnmapped Fastx --genomeDir /nfs/labs/nephew/human_databases/gencode_v25/ --readFilesIn $file --readFilesCommand zcat --outFilterMultimapNmax 25 --outFileNamePrefix $outdir/ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts";
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
  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 star.$f.sh");
  system("qsub -q bigmem -pe pe_slots 8 -wd $cwd star.$f.sh");
}
exit;

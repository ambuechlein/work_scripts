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
  next if($f =~ /Undetermined/o);
  next unless($f =~ /\.nonredundant.fasta$/o);
  my $file = $dir . '/' . $f;
  my $outdir = $out.'/'.$f.'_thout';

  my $command = 'tophat2 -o '.$outdir.' --transcriptome-index=/nfs/labs/nephew/human_databases/gencode_v25/transcriptome-index -p 8 --b2-very-sensitive --read-edit-dist 2 --max-multihits 100 --library-type fr-secondstrand /nfs/labs/nephew/human_databases/gencode_v25/genome ' . $file;
  my $command = "hisat2 -f --rna-strandness F -p 8 -x /nfs/labs/nephew/human_databases/gencode_v25/genome -U $file -S $outdir/accepted_hits.sam";
  open(my $sh, ">tophat.$f.sh") or die "can't open shell script $f: $!\n\n";

  print $sh '#!/bin/bash';
  print $sh "\n";
  print $sh 'echo ">>>>> Tophat "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'mkdir -p '.$outdir . ' || { echo "mkdir failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Tophat"', "\n";
  print $sh $command . ' || { echo "tophat failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Samtools"', "\n";
  print $sh '/nfs/bio/sw/bin/samtools sort -@ 8 -m 1G ' .$outdir. '/accepted_hits.bam ' . $outdir . '/accepted_hits_sorted' . ' || { echo "samtools sort failed"; exit 1; }'."\n";
  print $sh '/nfs/bio/sw/bin/samtools index '.$outdir. '/accepted_hits_sorted.bam' . ' || { echo "samtools index failed"; exit 1; }'."\n";
#  print $sh '/nfs/bio/sw/bin/samtools view ' .$outdir. '/accepted_hits_sorted.bam > ' . $outdir . '/accepted_hits_sorted.sam' . ' || { echo "samtools view failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 tophat.$f.sh");
  system("qsub -q bigmem,cluster -pe pe_slots 8 -wd $cwd tophat.$f.sh");
}
exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $results = GetOptions(\%options,
                         "directory|d=s",
                         "tophat|t=s",
                         "repeatmasker|r=s",
                         "summary|s=s",
                         "out|o=s",
                        );
my $pwd = cwd();
my $dir = $options{directory};

opendir(my $din, $dir) or die "can't open directory $dir: $!\n";
while(my $f = readdir($din)){
  next if($f =~ /Undetermined/o);
  next unless($f =~ /\.nonredundant.map$/o);
  my $file = $dir . '/' . $f;

#  $file =~ /(\w+)\.fastq\.nonredundant\.map$/go;
  $file =~ /((\w|-)+)\.(fastq(\.gz)?\.nonredundant)\.map$/go;
  my $col1 = $1;
  my $post = $3;
#  warn "$file\t$col1\t$post\n";
  open(my $dm, ">$col1.$options{out}") or die "can't open datamap output: $!\n";
  open(my $sh, ">combineresults.$col1.sh") or die "can't open shell script output: $!\n";

#  my $col2 = "$options{tophat}/$col1.fastq.nonredundant.fasta_thout/accepted_hits_sorted.bam";
  my $col2 = "$options{tophat}/$col1.$post.fasta_thout/accepted_hits_sorted.bam";
  die "$col2 doesn't exits" unless(-e $col2);
  my $col3 = "$dir/$f";
  die "$col3 doesn't exits" unless(-e $col3);
#  my $col4 = "$options{repeatmasker}/$col1.fastq.nonredundant.fasta.tsv";
  my $col4 = "$options{repeatmasker}/$col1.$post.fasta.tsv";
  die "$col4 doesn't exits" unless(-e $col4);
  my $col5 = "$options{summary}/$col1.counts.sorted.gz";
  die "$col5 doesn't exits" unless(-e $col5);
  my $col6 = "holder";
  my $col7 = "holder";
#  my $col8 = "$options{tophat}/$col1.fastq.nonredundant.fasta_thout/unmapped.bam";
  my $col8 = "$options{tophat}/$col1.$post.fasta_thout/unmapped.bam";
  die "$col8 doesn't exits" unless(-e $col8);
  my $col9 = "$options{summary}/$col1.sam.gz";
  my $col10 = "holder";
  my $col11 = "holder";

  print $dm "$col1\t$col2\t$col3\t$col4\t$col5\t$col6\t$col7\t$col8\t$col9\t$col10\t$col11\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> step3 combineDataSets "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n\n";
  print $sh "perl /nfs/labs/nephew/scripts/rnaseq_summary/combinedNewNephewDatasets2.fixTags.new.pl $pwd/$col1.$options{out}" . ' || { echo "combine failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $dm;
  close $sh;
  system("chmod 755 combineresults.$col1.sh");
  system('qsub -q bigmem -pe pe_slots 8 ' . "-wd $cwd combineresults.$col1.sh");
}
close $din;

exit;

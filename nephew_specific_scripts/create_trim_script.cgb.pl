#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my %options;
my $results = GetOptions(\%options,
                         "directory|d=s",
                         "map|m=s",
                         "outdirectory|o=s",
                        );
my $cwd = cwd();
my %lib = (
           'D701.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCTTCTGCTTG',
           'D702.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCGGAGAATCTCGTATGCCGTCTTCTGCTTG',
           'D703.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGCTCATTATCTCGTATGCCGTCTTCTGCTTG',
           'D704.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGATTCCATCTCGTATGCCGTCTTCTGCTTG',
           'D705.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCAGAAATCTCGTATGCCGTCTTCTGCTTG',
           'D706.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAATTCGTATCTCGTATGCCGTCTTCTGCTTG',
           'D707.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTGAAGCTATCTCGTATGCCGTCTTCTGCTTG',
           'D708.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAATGCGCATCTCGTATGCCGTCTTCTGCTTG',
           'D709.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGGCTATGATCTCGTATGCCGTCTTCTGCTTG',
           'D710.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCGCGAAATCTCGTATGCCGTCTTCTGCTTG',
           'D711.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCTCGCGCATCTCGTATGCCGTCTTCTGCTTG',
           'D712.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCGATAGATCTCGTATGCCGTCTTCTGCTTG',
           'D7XX.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATATCTCGTATGCCGTCTTCTGCTTG',
          );
my %bar2vec = (
               'ATTACTCG' => 'D701.vec.fa',
               'TCCGGAGA' => 'D702.vec.fa',
               'CGCTCATT' => 'D703.vec.fa',
               'GAGATTCC' => 'D704.vec.fa',
               'ATTCAGAA' => 'D705.vec.fa',
               'GAATTCGT' => 'D706.vec.fa',
               'CTGAAGCT' => 'D707.vec.fa',
               'TAATGCGC' => 'D708.vec.fa',
               'CGGCTATG' => 'D709.vec.fa',
               'TCCGCGAA' => 'D710.vec.fa',
               'TCTCGCGC' => 'D711.vec.fa',
               'AGCGATAG' => 'D712.vec.fa',
               'ATCACGAT' => 'D7XX.vec.fa',
);
open(my $min, $options{map}) or die "Can't open mapping file: $!\n";
<$min>;
my %sample2bar;
while(<$min>){
# Miller-19Nov2015-17smRNA_with_PEG	TTAGGC	22438216	34
  chomp $_;
  my ($l, $b) = split(/\t/, $_);
  $sample2bar{$l} = $b;
}
close $min;

my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
while(my $file = readdir($din)){
  next unless($file =~ /fastq.gz$/);
# GSF1021-RunA-H28-KTB84-NP-H _S28_R1_001.fastq.gz
  my $vec = '';
  if($file =~ /(\S+)_S\d+_R1_001\.fastq\.gz/){
    die "$file did not have an adapter\n" if(not defined $1);
#    $vec = "newvec$1.fa";
    $vec = $bar2vec{$sample2bar{$1}};
    warn "$file did not have an adapter file\n$1\t$sample2bar{$1} $bar2vec{$sample2bar{$1}}\n\n" if(not defined $bar2vec{$sample2bar{$1}});
  } else {
    warn "Houston we have a problem: $file\n";
  }

  my $prefix = $file;
  $prefix =~ s/_R1_001\.fastq\.gz$//o;
  
  my $newfile1 = "$options{outdirectory}/step1/$file";
  my $newfile2 = "$options{outdirectory}/step5/$file";
  my $newfile3 = "$options{outdirectory}/reverse_complement/$file";
#  my $newfile4 = "$options{outdirectory}/non_redundant/$file";

  $newfile1 =~ s/\.gz$//;
  $newfile2 =~ s/\.gz$//;
#  $newfile4 =~ s/\.gz$//;

  system("mkdir -p $options{outdirectory}/step1") if(not -e "$options{outdirectory}/step1");
  system("mkdir -p $options{outdirectory}/step5") if(not -e "$options{outdirectory}/step5");
  system("mkdir -p $options{outdirectory}/reverse_complement") if(not -e "$options{outdirectory}/reverse_complement");
#  system("mkdir -p $options{outdirectory}/non_redundant") if(not -e "$options{outdirectory}/non_redundant");

  my $command1 = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile1.' ILLUMINACLIP:/nfs/labs/nephew/adapters/d701_adapters/'.$vec.':2:20:5'.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17';

  my $command2 = "perl /home/abuechle/bin/truncateGReads.pl -i $newfile1 -o $newfile2 -l 18";
  my $command3 = "/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_reverse_complement -z -i $newfile2 -o $newfile3";
#  my $command4 = 'perl /nfs/labs/nephew/scripts/reduce_fastq_lib_sge.pl '.$newfile3. ' '.$newfile4;
  my $command5 = 'perl /home/abuechle/bin/countTrimmedReads.pl -a '.$dir.'/'.$file.' -b '.$newfile2.' -o '.$dir.'/'.$file.'.count -p '. $prefix;

  my $command6 = 'rm '.$newfile1;
  my $command7 = "pigz $newfile2";

  open(my $sh, ">trim.$file.fastx.sh") or die "Can't open shell script: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> Trimming "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting trimmomatic 1"', "\n";
  print $sh $command1 . ' || { echo "trimmomatic failed"; exit 1; }'."\n";
  print $sh 'echo "truncating Gs"', "\n";
  print $sh $command2 . ' || { echo "truncateG failed"; exit 1; }'."\n";
  print $sh 'echo "reverse complement"', "\n";
  print $sh $command3 . ' || { echo "reverse failed"; exit 1; }'."\n";
#  print $sh 'echo "reducing"', "\n";
#  print $sh $command4 . ' || { echo "reduce failed"; exit 1; }'."\n";
  print $sh 'echo "counting"', "\n";
  print $sh $command5 . ' || { echo "count1 failed"; exit 1; }'."\n";
  print $sh 'echo "gzip"', "\n";
  print $sh $command6 . ' || { echo "gzip2 failed"; exit 1; }'."\n";
  print $sh $command7 . ' || { echo "gzip9 failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trim.$file.fastx.sh");
  system("qsub -q bigmem -pe pe_slots 8 -wd $cwd trim.$file.fastx.sh");
}
closedir $din;
exit;

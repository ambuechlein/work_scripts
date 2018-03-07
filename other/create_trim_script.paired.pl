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
            'TruSeq_Adapter_Index_1.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_2.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_3.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_4.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_5.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_6.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_7.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_8.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_9.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_10.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_11.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_12.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_13.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_14.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_15.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_16.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_18.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_19.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_20.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_21.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_22.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_23.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_25.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_27.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG',
            'TruSeq_Adapter_Index_23.2.vec.fa' => 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGTTCTCGTATGCCGTCTTCTGCTTG',
          );
my %bar2vec = (
            'ATCACG' => 'TruSeq_Adapter_Index_1.vec.fa',
            'CGATGT' => 'TruSeq_Adapter_Index_2.vec.fa',
            'TTAGGC' => 'TruSeq_Adapter_Index_3.vec.fa',
            'TGACCA' => 'TruSeq_Adapter_Index_4.vec.fa',
            'ACAGTG' => 'TruSeq_Adapter_Index_5.vec.fa',
            'GCCAAT' => 'TruSeq_Adapter_Index_6.vec.fa',
            'CAGATC' => 'TruSeq_Adapter_Index_7.vec.fa',
            'ACTTGA' => 'TruSeq_Adapter_Index_8.vec.fa',
            'GATCAG' => 'TruSeq_Adapter_Index_9.vec.fa',
            'TAGCTT' => 'TruSeq_Adapter_Index_10.vec.fa',
            'GGCTAC' => 'TruSeq_Adapter_Index_11.vec.fa',
            'CTTGTA' => 'TruSeq_Adapter_Index_12.vec.fa',
            'AGTCAA' => 'TruSeq_Adapter_Index_13.vec.fa',
            'AGTTCC' => 'TruSeq_Adapter_Index_14.vec.fa',
            'ATGTCA' => 'TruSeq_Adapter_Index_15.vec.fa',
            'CCGTCC' => 'TruSeq_Adapter_Index_16.vec.fa',
            'GTCCGC' => 'TruSeq_Adapter_Index_18.vec.fa',
            'GTGAAA' => 'TruSeq_Adapter_Index_19.vec.fa',
            'GTGGCC' => 'TruSeq_Adapter_Index_20.vec.fa',
            'GTTTCG' => 'TruSeq_Adapter_Index_21.vec.fa',
            'CGTACG' => 'TruSeq_Adapter_Index_22.vec.fa',
            'CGTACG' => 'TruSeq_Adapter_Index_23.vec.fa',
            'ACTGAT' => 'TruSeq_Adapter_Index_25.vec.fa',
            'ATTCCT' => 'TruSeq_Adapter_Index_27.vec.fa',
            'GAGTGG' => 'TruSeq_Adapter_Index_23.2.vec.fa',
);
open(my $min, $options{map}) or die "Can't open mapping file: $!\n";
<$min>;
my %sample2bar;
while(<$min>){
  chomp $_;
# Sample	BarcodeSequence
# GSF1565A-BrainSBT-P10-11	AGTTCC
# GSF1565A-BrainSBT-P10-12	ATGTCA
# GSF1565A-BrainSBT-P10-14	GTGAAA
  my ($s, $b)  = split(/\t/, $_);
  $sample2bar{$s} = $b;
}
close $min;

my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
while(my $file = readdir($din)){
  next unless($file =~ /R1_001\.fastq\.gz$/);
# GSF1565A-BrainSBT-P10-6_S13_R1_001.fastq.gz
  my $vec = '';
  if($file =~ /(\S+)_S\d+_R1_001\.fastq\.gz/){
    my $s = $1;
    die "$file did not have an adapter\n" if(not defined $s);
    $vec = $bar2vec{$sample2bar{$s}};
    warn "$file ($s) did not have an adapter file\n$s\t$sample2bar{$s} $bar2vec{$sample2bar{$s}}\n\n" if(not defined $bar2vec{$sample2bar{$s}});

  } else {
    warn "Houston we have a problem: $file\n";
  }
##
  my $prefix = $file;
  $prefix =~ s/_L001_R1_001\.fastq\.gz$//o;
  my $file2 = $file;
  $file2 =~ s/R1_001\.fastq\.gz/R2_001\.fastq\.gz/go;
  my $newfile1 = "$options{outdirectory}/trimmed/$file";
  my $newfile2 = "$options{outdirectory}/trimmed/$file2";
  $newfile1 =~ s/\.fastq\.gz$//;
  $newfile2 =~ s/\.fastq\.gz$//;

  system("mkdir -p $options{outdirectory}/trimmed") if(not -e "$options{outdirectory}/trimmed");
  my $fiveprime = 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC';

  my $command1 = "java -jar /N/soft/rhel6/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 $dir/$file $dir/$file2 $newfile1.paired.fastq.gz $newfile1.unpaired.fastq.gz $newfile2.paired.fastq.gz $newfile2.unpaired.fastq.gz ILLUMINACLIP:/N/u/abuechle/Karst/bin/adapters/trueseq_adapters/$vec:2:20:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17";

  open(my $sh, ">trimpair.$file.sh") or die "Can't open shell script: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> Trimming "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting trimmomatic 1"', "\n";
  print $sh $command1 . ' || { echo "trimmomatic failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trimpair.$file.sh");
#  system("qsub -q bigmem -pe pe_slots 8 -wd $cwd trimpair.$file.sh");
}
closedir $din;
exit;

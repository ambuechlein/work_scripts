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
# Lane    Sample_ID       index   Sample_Project
# 6       152_UnTx_PDX    ATCACG  DaveMiller
# 6       146_4wk_MTD     CAGATC  DaveMiller
  chomp $_;
#  my ($l, $sp, $d, $i, $b) = split(/\t/, $_);
  my ($s, $b)  = split(/\t/, $_);
  $sample2bar{$s} = $b;
}
close $min;

my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
while(my $file = readdir($din)){
  next unless($file =~ /fastq.gz$/);
# 136-44_S56_L008_R1_001.fastq.gz
# 146-4wk-MTD_S34_L006_R1_001.fastq.gz
# 146-patient_S35_L006_R1_001.fastq.gz
# 152-4wk-MTD_S36_L006_R1_001.fastq.gz
# 152-UnTx-PDX_S33_L006_R1_001.fastq.gz
  my $vec = '';
  if($file =~ /(\S+)_S\d+_L\d+_R1_001\.fastq\.gz/){
    my $s = $1;
    $s =~ s/\-/_/go;
    die "$file did not have an adapter\n" if(not defined $s);
#    $vec = "newvec$1.fa";
    $vec = $bar2vec{$sample2bar{$s}};
    warn "$file did not have an adapter file\n$s\t$sample2bar{$s} $bar2vec{$sample2bar{$s}}\n\n" if(not defined $bar2vec{$sample2bar{$s}});
  } else {
    warn "Houston we have a problem: $file\n";
  }

  my $prefix = $file;
  $prefix =~ s/_R1_001\.fastq\.gz$//o;
  
#  my $newfile1 = "$options{outdirectory}/step1/$file";
  my $newfile2 = "$options{outdirectory}/step1/$file";
  my $newfile3 = "$options{outdirectory}/step2/$file";
  my $newfile4 = "$options{outdirectory}/step3/$file";
  my $newfile5 = "$options{outdirectory}/step4/$file";
  my $newfile6 = "$options{outdirectory}/non_redundant/$file";
  my $newfile7 = "$options{outdirectory}/step5/$file";
#  $newfile1 =~ s/\.gz$//;
  $newfile2 =~ s/\.gz$//;
  $newfile3 =~ s/\.gz$//;
  $newfile4 =~ s/\.gz$//;
  $newfile5 =~ s/\.gz$//;
  $newfile6 =~ s/\.gz$//;
  $newfile7 =~ s/\.gz$//;
  my $newfile3_trimmed = "$newfile3.trimmed.fastq";
  my $newfile3_untrimmed = "$newfile3.untrimmed.fastq";

  system("mkdir -p $options{outdirectory}/step1") if(not -e "$options{outdirectory}/step1");
  system("mkdir -p $options{outdirectory}/step2") if(not -e "$options{outdirectory}/step2");
  system("mkdir -p $options{outdirectory}/step3") if(not -e "$options{outdirectory}/step3");
  system("mkdir -p $options{outdirectory}/step4") if(not -e "$options{outdirectory}/step4");
  system("mkdir -p $options{outdirectory}/step5") if(not -e "$options{outdirectory}/step5");
#  system("mkdir -p $options{outdirectory}/step6") if(not -e "$options{outdirectory}/step6");
  system("mkdir -p $options{outdirectory}/non_redundant") if(not -e "$options{outdirectory}/non_redundant");


#  my $command = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile1.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17';
#  my $command1 = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$newfile1.' '.$newfile2.' ILLUMINACLIP:/nfs/labs/nephew/adapters/trueseq_adapters/'.$vec.':2:20:5';

  my $command = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile2.' ILLUMINACLIP:/nfs/labs/nephew/adapters/trueseq_adapters/'.$vec.':2:20:5'.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17';

  my $command2 = "perl /home/abuechle/bin/split_untrimmed.pl -a $dir/$file -b $newfile2 -o $newfile3";
  my $command3 = "perl /home/abuechle/bin/trimShortPotentialAdapters.pl $newfile3_untrimmed $lib{$vec} 10 2 17 > $newfile3.tmp";
  my $command4 = "cat $newfile3.tmp $newfile3_trimmed > $newfile3";

  my $fiveprime = 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC';
  my $command5 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$fiveprime.' -l 17 -v -M 1 -i '.$newfile3.' -o '.$newfile4;

  my $command6 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$lib{$vec}.' -l 17 -v -M 1 -i '.$newfile4.' -o '.$newfile5;
  my $command18 = "perl /home/abuechle/bin/truncateGReads.pl -i $newfile5 -o $newfile7 -l 18";
  my $command7 = 'perl /nfs/labs/nephew/all_new_data/fastq/non_redundant/reduce_fastq_lib_sge.pl '.$newfile7. ' '.$newfile6;
  my $command8 = 'perl /home/abuechle/bin/countTrimmedReads.pl -a '.$dir.'/'.$file.' -b '.$newfile7.' -o '.$dir.'/'.$file.'.count -p '. $prefix;

  my $command11 = 'pigz '.$newfile2;
  my $command12 = 'pigz '.$newfile3;
  my $command13 = 'pigz '.$newfile4;
  my $command14 = 'pigz '.$newfile5;
  my $command15 = 'pigz '.$newfile3_trimmed;
  my $command16 = 'pigz '.$newfile3_untrimmed;
  my $command17 = 'pigz '.$newfile3.'.tmp';
  my $command19 = 'pigz '.$newfile7;

  open(my $sh, ">trim.$file.fastx.sh") or die "Can't open shell script: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> Trimming "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting trimmomatic 1"', "\n";
  print $sh $command . ' || { echo "trimmomatic failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting trimmomatic 2"', "\n";
#  print $sh $command1 . ' || { echo "trimmomatic 2 failed"; exit 1; }'."\n";
  print $sh 'echo "splitting untrimmed"', "\n";
  print $sh $command2 . ' || { echo "split failed failed"; exit 1; }'."\n";
  print $sh 'echo "doug trimming"', "\n";
  print $sh $command3 . ' || { echo "doug trim failed"; exit 1; }'."\n";
  print $sh $command4 . ' || { echo "cat failed"; exit 1; }'."\n";
  print $sh 'echo "fastx1 trimming"', "\n";
  print $sh $command5 . ' || { echo "fastx1 failed"; exit 1; }'."\n";
  print $sh 'echo "fastx2 trimming"', "\n";
  print $sh $command6 . ' || { echo "fastx2 failed"; exit 1; }'."\n";
  print $sh 'echo "truncating Gs"', "\n";
  print $sh $command18 . ' || { echo "truncateG failed"; exit 1; }'."\n";
  print $sh 'echo "reducing"', "\n";
  print $sh $command7 . ' || { echo "reduce failed"; exit 1; }'."\n";
  print $sh 'echo "counting"', "\n";
  print $sh $command8 . ' || { echo "count1 failed"; exit 1; }'."\n";
  print $sh 'echo "gzip"', "\n";
#  print $sh $command10 . ' || { echo "gzip1 failed"; exit 1; }'."\n";
  print $sh $command11 . ' || { echo "gzip2 failed"; exit 1; }'."\n";
  print $sh $command12 . ' || { echo "gzip3 failed"; exit 1; }'."\n";
  print $sh $command13 . ' || { echo "gzip4 failed"; exit 1; }'."\n";
  print $sh $command14 . ' || { echo "gzip5 failed"; exit 1; }'."\n";
  print $sh $command15 . ' || { echo "gzip6 failed"; exit 1; }'."\n";
  print $sh $command16 . ' || { echo "gzip7 failed"; exit 1; }'."\n";
  print $sh $command17 . ' || { echo "gzip8 failed"; exit 1; }'."\n";
  print $sh $command19 . ' || { echo "gzip9 failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trim.$file.fastx.sh");
  system("qsub -q bigmem\@intron,bigmem\@exon.cgb.indiana.edu -pe pe_slots 4 -wd $cwd trim.$file.fastx.sh");
}
closedir $din;
exit;

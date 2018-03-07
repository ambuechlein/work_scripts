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
                         "reduce|r=s",
                        );
if(not defined $options{reduce}){
  $options{reduce} = 1;
}
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
# Sample_Name	I7_Index_ID	index
# GSF1044A-D1-Flag-2	A002	ATCACG
# GSF1044A-D2-Flag-3	A003	TTAGGC
  chomp $_;
#  my ($l, $sp, $d, $i, $b) = split(/\t/, $_);
  my ($s, $b)  = split(/\t/, $_);
  $sample2bar{$s} = $b;
}
close $min;

my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
while(my $file = readdir($din)){
#  next unless($file =~ /^GSF1348A\-redo/);
  next unless($file =~ /fastq.gz$/);
# GSF1118-MCF-7-Hypo-p2_S2_R1_001
  my $vec = '';
  if($file =~ /(\S+)_S\d+_R1_001\.fastq\.gz/){
    my $s = $1;
#    $s =~ s/\-/_/go;
    die "$file did not have an adapter\n" if(not defined $s);
    $vec = $bar2vec{$sample2bar{$s}};
    warn "$file did not have an adapter file\n$s\t$sample2bar{$s} $bar2vec{$sample2bar{$s}}\n\n" if(not defined $bar2vec{$sample2bar{$s}});
  } else {
    warn "Houston we have a problem: $file\n";
  }

  my $prefix = $file;
  $prefix =~ s/_R1_001\.fastq\.gz$//o;
  
  my $newfile1 = "$options{outdirectory}/step1/$file";
  my $newfile2 = "$options{outdirectory}/step2/$file";
  my $newfile3 = "$options{outdirectory}/step3/$file";
  my $newfile4 = "$options{outdirectory}/step4/$file";
  my $newfile5 = "$options{outdirectory}/step5/$file";
  my $newfile6 = "$options{outdirectory}/reverse_complement/$file";
  my $newfile7 = "$options{outdirectory}/non_redundant/$file";

  $newfile1 =~ s/\.gz$//;
  $newfile2 =~ s/\.gz$//;
  $newfile3 =~ s/\.gz$//;
  $newfile4 =~ s/\.gz$//;
  $newfile5 =~ s/\.gz$//;
#  $newfile6 =~ s/\.gz$//;
  my $newfile2_trimmed = "$newfile2.trimmed.fastq";
  my $newfile2_untrimmed = "$newfile2.untrimmed.fastq";

  system("mkdir -p $options{outdirectory}/step1") if(not -e "$options{outdirectory}/step1");
  system("mkdir -p $options{outdirectory}/step2") if(not -e "$options{outdirectory}/step2");
  system("mkdir -p $options{outdirectory}/step3") if(not -e "$options{outdirectory}/step3");
  system("mkdir -p $options{outdirectory}/step4") if(not -e "$options{outdirectory}/step4");
  system("mkdir -p $options{outdirectory}/step5") if(not -e "$options{outdirectory}/step5");
  system("mkdir -p $options{outdirectory}/reverse_complement") if(not -e "$options{outdirectory}/reverse_complement");
  system("mkdir -p $options{outdirectory}/non_redundant") if(not -e "$options{outdirectory}/non_redundant");

  my $fiveprime = 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC';

  my $command1 = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile1.' ILLUMINACLIP:/nfs/labs/nephew/adapters/trueseq_adapters/'.$vec.':2:20:5'.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17';

  my $command2 = "perl /home/abuechle/bin/split_untrimmed.pl -a $dir/$file -b $newfile1 -o $newfile2";
  my $command3 = "perl /home/abuechle/bin/trimShortPotentialAdapters.pl $newfile2_untrimmed $lib{$vec} 10 2 17 > $newfile2.tmp";
  my $command4 = "cat $newfile2.tmp $newfile2_trimmed > $newfile2";

  my $command5 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.14/bin/fastx_clipper -Q33 -a '.$fiveprime.' -l 17 -v -M 1 -i '.$newfile2.' -o '.$newfile3;

  my $command6 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.14/bin/fastx_clipper -Q33 -a '.$lib{$vec}.' -l 17 -v -M 1 -i '.$newfile3.' -o '.$newfile4;
  my $command7 = "perl /home/abuechle/bin/truncateGReads.pl -i $newfile4 -o $newfile5 -l 18";

  my $command8 = "/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_reverse_complement -Q33 -z -i $newfile5 -o $newfile6";

  my $command9 = 'perl /nfs/labs/nephew/scripts/reduce_fastq_lib_sge.pl '.$newfile6. ' '.$newfile7;
  my $command10 = 'perl /home/abuechle/bin/countTrimmedReads.pl -a '.$dir.'/'.$file.' -b '.$newfile5.' -o '.$dir.'/'.$file.'.count -p '. $prefix;

  my $command11 = 'pigz '.$newfile1;
  my $command12 = 'pigz '.$newfile2;
  my $command13 = 'pigz '.$newfile2_trimmed;
  my $command14 = 'pigz '.$newfile2_untrimmed;
  my $command15 = 'pigz '.$newfile2.'.tmp';
  my $command16 = 'pigz '.$newfile3;
  my $command17 = 'pigz '.$newfile4;
  my $command18 = 'pigz '.$newfile5;

  open(my $sh, ">trim.$file.sh") or die "Can't open shell script: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> Trimming "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting trimmomatic 1"', "\n";
  print $sh $command1 . ' || { echo "trimmomatic failed"; exit 1; }'."\n";
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
  print $sh $command7 . ' || { echo "truncateG failed"; exit 1; }'."\n";
  print $sh 'echo "reverse complement"', "\n" if($options{reduce});
  print $sh $command8 . ' || { echo "reverse complement failed"; exit 1; }'."\n" if($options{reduce});
  print $sh 'echo "reducing"', "\n" if($options{reduce});
  print $sh $command9 . ' || { echo "reduce failed"; exit 1; }'."\n" if($options{reduce});
  print $sh 'echo "counting"', "\n";
  print $sh $command10 . ' || { echo "counting failed"; exit 1; }'."\n";
  print $sh 'echo "gzip"', "\n";
  print $sh $command11 . ' || { echo "gzip1 failed"; exit 1; }'."\n";
  print $sh $command12 . ' || { echo "gzip2 failed"; exit 1; }'."\n";
  print $sh $command13 . ' || { echo "gzip3 failed"; exit 1; }'."\n";
  print $sh $command14 . ' || { echo "gzip4 failed"; exit 1; }'."\n";
  print $sh $command15 . ' || { echo "gzip5 failed"; exit 1; }'."\n";
  print $sh $command16 . ' || { echo "gzip6 failed"; exit 1; }'."\n";
  print $sh $command17 . ' || { echo "gzip7 failed"; exit 1; }'."\n";
  print $sh $command18 . ' || { echo "gzip8 failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trim.$file.sh");
#  system("qsub -q bigmem\@exon.cgb.indiana.edu,bigmem\@intron,bigmem\@antibody,bigmem\@antigen -pe pe_slots 4 -wd $cwd trim.$file.sh");
}
closedir $din;
exit;

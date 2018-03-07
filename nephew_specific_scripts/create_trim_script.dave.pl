#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions(\%options,
                         "directory|d=s",
                         "map|m=s",
                         "outdirectory|o=s",
                        );
my %lib = (
            'newvec.fa'  => 'CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec1.fa' => 'CAAGCAGAAGACGGCATACGAGATATCACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec2.fa' => 'CAAGCAGAAGACGGCATACGAGATCGATGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec3.fa' => 'CAAGCAGAAGACGGCATACGAGATTTAGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec4.fa' => 'CAAGCAGAAGACGGCATACGAGATTGACCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec5.fa' => 'CAAGCAGAAGACGGCATACGAGATACAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec6.fa' => 'CAAGCAGAAGACGGCATACGAGATGCCAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec7.fa' => 'CAAGCAGAAGACGGCATACGAGATCAGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
            'newvec8.fa' => 'CAAGCAGAAGACGGCATACGAGATACTTGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
          );
my %bar2vec = (
           'ATCACG' => 'newvec1.fa',
           'CGATGT' => 'newvec2.fa',
           'TTAGGC' => 'newvec3.fa',
           'TGACCA' => 'newvec4.fa',
           'ACAGTG' => 'newvec5.fa',
           'GCCAAT' => 'newvec6.fa',
           'CAGATC' => 'newvec7.fa',
           'ACTTGA' => 'newvec8.fa',
);
open(my $min, $options{map}) or die "Can't open mapping file: $!\n";
<$min>;
my %sample2bar;
while(<$min>){
# Sample	Content	Qbit	Index	Index Seq
# 1	Wt-EZH2	50.4	1	ATCACG
# 2	Wt-EZH2	51.8	2	CGATGT
# 3	EZH2-T372E	46.2	3	TTAGGC

  chomp $_;
  my ($sp, $c, $q, $i, $b) = split(/\t/, $_);
  $sample2bar{$c} = $b;
}
close $min;

my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
while(my $file = readdir($din)){
  next unless($file =~ /fastq.gz$/);
# GSF844-1-Wt-EZH2_S1_R1_001.fastq.gz

  my $vec = 'newvec.fa';
  if($file =~ /^GSF844-(\d+)-\S+_R1_001\.fastq\.gz/){
    die "$file did not have an adapter\n" if(not defined $1);
    $vec = "newvec$1.fa";
#    $vec = $bar2vec{$sample2bar{$1}};
#    warn "$file did not have an adapter file\n$1\t$sample2bar{$1}\n\n" if(not defined $bar2vec{$sample2bar{$1}});
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
  my $newfile6 = "$options{outdirectory}/non_redundant/$file"; 
  $newfile1 =~ s/\.gz$//;
  $newfile2 =~ s/\.gz$//;
  $newfile3 =~ s/\.gz$//;
  $newfile4 =~ s/\.gz$//;
  $newfile5 =~ s/\.gz$//;
  $newfile6 =~ s/\.gz$//;
  my $newfile3_trimmed = "$newfile3.trimmed.fastq";
  my $newfile3_untrimmed = "$newfile3.untrimmed.fastq";

  system("mkdir -p $options{outdirectory}/step1") if(not -e "$options{outdirectory}/step1");
  system("mkdir -p $options{outdirectory}/step2") if(not -e "$options{outdirectory}/step2");
  system("mkdir -p $options{outdirectory}/step3") if(not -e "$options{outdirectory}/step3");
  system("mkdir -p $options{outdirectory}/step4") if(not -e "$options{outdirectory}/step4");
  system("mkdir -p $options{outdirectory}/step5") if(not -e "$options{outdirectory}/step5");
  system("mkdir -p $options{outdirectory}/non_redundant") if(not -e "$options{outdirectory}/non_redundant");


  my $command = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile1.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17';
  my $command1 = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$newfile1.' '.$newfile2.' ILLUMINACLIP:/nfs/labs/nephew/adapters/'.$vec.':2:20:5';

  my $command2 = "perl /home/abuechle/bin/split_untrimmed.pl -a $newfile1 -b $newfile2 -o $newfile3";
  my $command3 = "perl /home/abuechle/bin/trimShortPotentialAdapters.pl $newfile3_untrimmed $lib{$vec} 10 2 17 > $newfile3.tmp";
  my $command4 = "cat $newfile3.tmp $newfile3_trimmed > $newfile3";

  my $fiveprime = 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC';
  my $command5 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$fiveprime.' -l 17 -v -M 1 -i '.$newfile3.' -o '.$newfile4;

  my $command6 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$lib{$vec}.' -l 17 -v -M 1 -i '.$newfile4.' -o '.$newfile5;
  my $command7 = 'perl /nfs/labs/nephew/all_new_data/fastq/non_redundant/reduce_fastq_lib_sge.pl '.$newfile5. ' '.$newfile6;
  my $command8 = 'perl /home/abuechle/bin/countTrimmedReads.pl -a '.$dir.'/'.$file.' -b '.$newfile5.' -o '.$dir.'/'.$file.'.count -p '. $prefix;

  my $command10 = 'pigz '.$newfile1;
  my $command11 = 'pigz '.$newfile2;
  my $command12 = 'pigz '.$newfile3;
  my $command13 = 'pigz '.$newfile4;
  my $command14 = 'pigz '.$newfile5;
  my $command15 = 'pigz '.$newfile3_trimmed;
  my $command16 = 'pigz '.$newfile3_untrimmed;
  my $command17 = 'pigz '.$newfile3.'.tmp';

  open(my $sh, ">trim.$file.fastx.sh") or die "Can't open shell script: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh 'echo ">>>>> Trimming "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting trimmomatic 1"', "\n";
  print $sh $command . ' || { echo "trimmomatic failed"; exit 1; }'."\n";
  print $sh 'echo "Starting trimmomatic 2"', "\n";
  print $sh $command1 . ' || { echo "trimmomatic 2 failed"; exit 1; }'."\n";
  print $sh 'echo "splitting untrimmed"', "\n";
  print $sh $command2 . ' || { echo "split failed failed"; exit 1; }'."\n";
  print $sh 'echo "doug trimming"', "\n";
  print $sh $command3 . ' || { echo "doug trim failed"; exit 1; }'."\n";
  print $sh $command4 . ' || { echo "cat failed"; exit 1; }'."\n";
  print $sh 'echo "fastx1 trimming"', "\n";
  print $sh $command5 . ' || { echo "fastx1 failed"; exit 1; }'."\n";
  print $sh 'echo "fastx2 trimming"', "\n";
  print $sh $command6 . ' || { echo "fastx2 failed"; exit 1; }'."\n";
  print $sh 'echo "reducing"', "\n";
  print $sh $command7 . ' || { echo "reduce failed"; exit 1; }'."\n";
  print $sh 'echo "counting"', "\n";
  print $sh $command8 . ' || { echo "count1 failed"; exit 1; }'."\n";
  print $sh 'echo "gzip"', "\n";
  print $sh $command10 . ' || { echo "gzip1 failed"; exit 1; }'."\n";
  print $sh $command11 . ' || { echo "gzip2 failed"; exit 1; }'."\n";
  print $sh $command12 . ' || { echo "gzip3 failed"; exit 1; }'."\n";
  print $sh $command13 . ' || { echo "gzip4 failed"; exit 1; }'."\n";
  print $sh $command14 . ' || { echo "gzip5 failed"; exit 1; }'."\n";
  print $sh $command15 . ' || { echo "gzip6 failed"; exit 1; }'."\n";
  print $sh $command16 . ' || { echo "gzip7 failed"; exit 1; }'."\n";
  print $sh $command17 . ' || { echo "gzip8 failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trim.$file.fastx.sh");
  system('qsub -q bigmem -pe pe_slots 16 '."trim.$file.fastx.sh");
}
closedir $din;
exit;

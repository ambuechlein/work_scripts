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
           'RNAPCRPrimerIndex1.fsa'  => 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex2.fsa'  => 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex3.fsa'  => 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex4.fsa'  => 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex5.fsa'  => 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex6.fsa'  => 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex7.fsa'  => 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex8.fsa'  => 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex9.fsa'  => 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex10.fsa' => 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex11.fsa' => 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex12.fsa' => 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex13.fsa' => 'CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex14.fsa' => 'CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex15.fsa' => 'CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex16.fsa' => 'CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex17.fsa' => 'CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex18.fsa' => 'CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex19.fsa' => 'CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex20.fsa' => 'CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex21.fsa' => 'CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex22.fsa' => 'CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex23.fsa' => 'CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex24.fsa' => 'CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex25.fsa' => 'CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex26.fsa' => 'CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex27.fsa' => 'CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex28.fsa' => 'CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex29.fsa' => 'CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex30.fsa' => 'CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex31.fsa' => 'CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex32.fsa' => 'CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex33.fsa' => 'CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex34.fsa' => 'CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex35.fsa' => 'CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex36.fsa' => 'CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex37.fsa' => 'CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex38.fsa' => 'CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex39.fsa' => 'CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex40.fsa' => 'CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex41.fsa' => 'CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex42.fsa' => 'CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex43.fsa' => 'CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex44.fsa' => 'CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex45.fsa' => 'CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex46.fsa' => 'CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex47.fsa' => 'CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
           'RNAPCRPrimerIndex48.fsa' => 'CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA',
          );
my %bar2vec = (
           'ATCACG' => 'RNAPCRPrimerIndex1.fsa',
           'CGATGT' => 'RNAPCRPrimerIndex2.fsa',
           'TTAGGC' => 'RNAPCRPrimerIndex3.fsa',
           'TGACCA' => 'RNAPCRPrimerIndex4.fsa',
           'ACAGTG' => 'RNAPCRPrimerIndex5.fsa',
           'GCCAAT' => 'RNAPCRPrimerIndex6.fsa',
           'CAGATC' => 'RNAPCRPrimerIndex7.fsa',
           'ACTTGA' => 'RNAPCRPrimerIndex8.fsa',
           'GATCAG' => 'RNAPCRPrimerIndex9.fsa',
           'TAGCTT' => 'RNAPCRPrimerIndex10.fsa',
           'GGCTAC' => 'RNAPCRPrimerIndex11.fsa',
           'CTTGTA' => 'RNAPCRPrimerIndex12.fsa',
           'AGTCAA' => 'RNAPCRPrimerIndex13.fsa',
           'AGTTCC' => 'RNAPCRPrimerIndex14.fsa',
           'ATGTCA' => 'RNAPCRPrimerIndex15.fsa',
           'CCGTCC' => 'RNAPCRPrimerIndex16.fsa',
           'GTAGAG' => 'RNAPCRPrimerIndex17.fsa',
           'GTCCGC' => 'RNAPCRPrimerIndex18.fsa',
           'GTGAAA' => 'RNAPCRPrimerIndex19.fsa',
           'GTGGCC' => 'RNAPCRPrimerIndex20.fsa',
           'GTTTCG' => 'RNAPCRPrimerIndex21.fsa',
           'CGTACG' => 'RNAPCRPrimerIndex22.fsa',
           'GAGTGG' => 'RNAPCRPrimerIndex23.fsa',
           'GGTAGC' => 'RNAPCRPrimerIndex24.fsa',
           'ACTGAT' => 'RNAPCRPrimerIndex25.fsa',
           'ATGAGC' => 'RNAPCRPrimerIndex26.fsa',
           'ATTCCT' => 'RNAPCRPrimerIndex27.fsa',
           'CAAAAG' => 'RNAPCRPrimerIndex28.fsa',
           'CAACTA' => 'RNAPCRPrimerIndex29.fsa',
           'CACCGG' => 'RNAPCRPrimerIndex30.fsa',
           'CACGAT' => 'RNAPCRPrimerIndex31.fsa',
           'CACTCA' => 'RNAPCRPrimerIndex32.fsa',
           'CAGGCG' => 'RNAPCRPrimerIndex33.fsa',
           'CATGGC' => 'RNAPCRPrimerIndex34.fsa',
           'CATTTT' => 'RNAPCRPrimerIndex35.fsa',
           'CCAACA' => 'RNAPCRPrimerIndex36.fsa',
           'CGGAAT' => 'RNAPCRPrimerIndex37.fsa',
           'CTAGCT' => 'RNAPCRPrimerIndex38.fsa',
           'CTATAC' => 'RNAPCRPrimerIndex39.fsa',
           'CTCAGA' => 'RNAPCRPrimerIndex40.fsa',
           'GACGAC' => 'RNAPCRPrimerIndex41.fsa',
           'TAATCG' => 'RNAPCRPrimerIndex42.fsa',
           'TACAGC' => 'RNAPCRPrimerIndex43.fsa',
           'TATAAT' => 'RNAPCRPrimerIndex44.fsa',
           'TCATTC' => 'RNAPCRPrimerIndex45.fsa',
           'TCCCGA' => 'RNAPCRPrimerIndex46.fsa',
           'TCGAAG' => 'RNAPCRPrimerIndex47.fsa',
           'TCGGCA' => 'RNAPCRPrimerIndex48.fsa',
);
open(my $min, $options{map}) or die "Can't open mapping file: $!\n";
<$min>;
my %sample2bar;
while(<$min>){
# Name on tube	SmallRNA index code
# 2746	ACTGAT
# 3697	ATGAGC
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
  next unless($file =~ /small/go);
# GSF1016-smallRNA-2717_S25_R1_001.fastq.gz
# GSF1016-smallRNA-2746_S19_R1_001.fastq.gz
  my $vec = '';
  if($file =~ /(GSF\S+)_S\d+_R1_001\.fastq\.gz/){
    my $s = $1;
#    $s =~ s/\-/_/go;
    die "$file did not have an adapter\n" if(not defined $s);
    $vec = $bar2vec{$sample2bar{$s}};
#    warn "$file did not have an adapter file\n$s\t$sample2bar{$s} $bar2vec{$sample2bar{$s}}\n\n" if(not defined $bar2vec{$sample2bar{$s}});
    warn "$s\t$sample2bar{$s}\n" if(not defined $bar2vec{$sample2bar{$s}});
next if(not defined $bar2vec{$sample2bar{$s}});
  } else {
    warn "Houston we have a problem: $file\n";
next;
  }
die "Adapter not found $vec\n" unless(-e "/nfs/labs/nephew/adapters/small_illumina/$vec");
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

  my $fiveprime = 'GTTCAGAGTTCTACAGTCCGACGATC';

  my $command1 = 'java -classpath /home/abuechle/bin/Trimmomatic-0.32_patched/dist/jar/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 16 -phred33 '.$dir.'/'.$file.' '.$newfile1.' ILLUMINACLIP:/nfs/labs/nephew/adapters/small_illumina/'.$vec.':2:20:5'.' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16';

  my $command2 = "perl /home/abuechle/bin/split_untrimmed.pl -a $dir/$file -b $newfile1 -o $newfile2";
  my $command3 = "perl /home/abuechle/bin/trimShortPotentialAdapters.pl $newfile2_untrimmed $lib{$vec} 10 2 16 > $newfile2.tmp";
  my $command4 = "cat $newfile2.tmp $newfile2_trimmed > $newfile2";

  my $command5 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$fiveprime.' -l 16 -v -M 1 -i '.$newfile2.' -o '.$newfile3;

  my $command6 = '/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_clipper -a '.$lib{$vec}.' -l 16 -v -M 1 -i '.$newfile3.' -o '.$newfile4;
  my $command7 = "perl /home/abuechle/bin/truncateGReads.pl -i $newfile4 -o $newfile5 -l 16";

  my $command8 = "/home/abuechle/bin/fastx_toolkit/fastx_toolkit-0.0.13.2/bin/fastx_reverse_complement -z -i $newfile5 -o $newfile6";

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
  print $sh 'echo "reverse complement"', "\n";
  print $sh $command8 . ' || { echo "reduce failed"; exit 1; }'."\n";
  print $sh 'echo "reducing"', "\n";
  print $sh $command9 . ' || { echo "count1 failed"; exit 1; }'."\n";
  print $sh 'echo "counting"', "\n";
  print $sh $command10 . ' || { echo "gzip1 failed"; exit 1; }'."\n";
  print $sh 'echo "gzip"', "\n";
  print $sh $command11 . ' || { echo "gzip2 failed"; exit 1; }'."\n";
  print $sh $command12 . ' || { echo "gzip3 failed"; exit 1; }'."\n";
  print $sh $command13 . ' || { echo "gzip4 failed"; exit 1; }'."\n";
  print $sh $command14 . ' || { echo "gzip5 failed"; exit 1; }'."\n";
  print $sh $command15 . ' || { echo "gzip6 failed"; exit 1; }'."\n";
  print $sh $command16 . ' || { echo "gzip7 failed"; exit 1; }'."\n";
  print $sh $command17 . ' || { echo "gzip8 failed"; exit 1; }'."\n";
  print $sh $command18 . ' || { echo "gzip9 failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";

  close $sh;
  system("chmod 755 trim.$file.sh");
  system("qsub -q bigmem -pe pe_slots 4 -wd $cwd trim.$file.sh");
}
closedir $din;
exit;

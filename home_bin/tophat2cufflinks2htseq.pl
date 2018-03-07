#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions (\%options,
                          "directory|d=s",
                          "output|o=s",
                         );
my $dir = $options{directory};
my $outdir = $options{output};

opendir(my $din, $dir) or die "can't open directroy: $!\n";
while(my $f = readdir($din)){
  next unless($f =~ /\.fastq$/o);
  my $file = $dir . '/' . $f;

  my $command = "tophat2 -o ${outdir}/${f}_thout --transcriptome-index=/nfs/projects/solexa/analysis/human.pomerening/tophat/transcriptome_index/big -p 8 --b2-very-sensitive --read-edit-dist 2 --max-multihits 100 /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/genome $file";

  open(my $sh, ">tophat.$f.sh") or die "can't open shell script $f: $!\n\n";

  print $sh '#!/bin/bash';
  print $sh "\n";
  print $sh 'echo ">>>>> readmapping with Tophat "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh "mkdir -p $outdir/${f}_thout/" . ' || { echo "mkdir failed"; exit 1; }'."\n";
  print $sh 'echo "Starting Tophat"', "\n";
  print $sh $command . ' || { echo "tophat failed"; exit 1; }'."\n";
  print $sh "/nfs/bio/sw/bin/samtools sort  $outdir/${f}_thout/accepted_hits.bam $outdir/${f}_thout/accepted_hits_sorted" . ' || { echo "samtools sort failed"; exit 1; }'."\n";
  print $sh "/nfs/bio/sw/bin/samtools index $outdir/${f}_thout/accepted_hits_sorted.bam" . ' || { echo "samtools index failed"; exit 1; }'."\n";
  print $sh "/nfs/bio/sw/bin/samtools view  $outdir/${f}_thout/accepted_hits_sorted.bam > $outdir/${f}_thout/accepted_hits_sorted.sam" . ' || { echo "samtools view failed"; exit 1; }'."\n";

  print $sh 'echo "Starting Stats"', "\n";
  print $sh "/nfs/bio/sw/bin/samtools flagstat $outdir/${f}_thout/accepted_hits_sorted.bam > $outdir/${f}_thout/accepted_hits_sorted.bam.stats" . ' || { echo "samtools flagstat failed"; exit 1; }'."\n";
#  print $sh 'READ1=`wc -l '.$file.' | gawk \'{print int($1/4)}\' `' . ' || { echo "gawk failed"; exit 1; }'."\n";
  print $sh 'READ1=`wc -l '.$file.' | awk \'{print int($1/4)}\' `' . ' || { echo "awk failed"; exit 1; }'."\n";
  print $sh 'FASTQREADS=$READ1', "\n";
  print $sh 'echo $FASTQREADS" fastq reads" >> '."$outdir/${f}_thout/accepted_hits_sorted.bam.stats" . ' || { echo "echo reads failed"; exit 1; }'."\n";
  print $sh 'JUNCTION=$(wc -l '.$outdir.'/'.$f.'_thout/junctions.bed | cut -d\' \' -f 1)' . ' || { echo "junction count failed"; exit 1; }'."\n";
  print $sh 'echo $JUNCTION" junction reads" >> '."$outdir/${f}_thout/accepted_hits_sorted.bam.stats" . ' || { echo "echo junctions failed"; exit 1; }'."\n";
  print $sh 'JUNCTGENE=$(/home/abuechle/bin/bedtools-2.17.0/bin/windowBed -a '.$outdir.'/'.$f.'_thout/junctions.bed -b /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/gencode.v16.annotation.gtf -u -w 200 | wc -l | cut -d\' \' -f 1)' . ' || { echo "echo junction windowBed failed"; exit 1; }'."\n";
  print $sh 'echo $JUNCTGENE" junction reads Gencode" >> '."$outdir/${f}_thout/accepted_hits_sorted.bam.stats" . ' || { echo "echo junctions2 failed"; exit 1; }'."\n";
#  print $sh 'echo "Starting Samstat"', "\n";
#  print $sh "samstat $outdir/${f}_thout/accepted_hits_sorted.bam" . ' || { echo "samstat failed"; exit 1; }'."\n";

  print $sh 'echo "Starting Cufflinks"', "\n";
  print $sh "mkdir -p $outdir/${f}_clout/" . ' || { echo "mkdir failed"; exit 1; }'."\n";
#if($useGencode){
  print $sh "cufflinks --quiet --GTF-guide /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/gencode.v16.annotation.gtf -p 8 -o $outdir/${f}_clout/ $outdir/${f}_thout/accepted_hits_sorted.bam" . ' || { echo "cufflinks failed"; exit 1; }'."\n";
#}else{
    # non reference guided
#  print $sh "cufflinks --quiet --frag-bias-correct /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/genome.fa -p 8 -o $outdir/${f}_clout/ $outdir/${f}_thout/accepted_hits_sorted.bam" . ' || { echo "cufflinks failed"; exit 1; }'."\n";
#}

  print $sh 'echo "Starting Read Counting"', "\n";
  print $sh "htseq-count --quiet $outdir/${f}_thout/accepted_hits_sorted.sam $outdir/${f}_clout/transcripts.gtf > $outdir/${f}_thout/cufflinks.gene" . ' || { echo "htseq cufflinks failed"; exit 1; }'."\n";
  print $sh "htseq-count --quiet --idattr=\"transcript_id\" $outdir/${f}_thout/accepted_hits_sorted.sam $outdir/${f}_clout/transcripts.gtf > $outdir/${f}_thout/cufflinks.transcripts" . ' || { echo "htseq cufflinks2 failed"; exit 1; }'."\n";

  print $sh "htseq-count --quiet $outdir/${f}_thout/accepted_hits_sorted.sam /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/gencode.v16.annotation.gtf | grep ENSG > $outdir/${f}_thout/gencode.gene" . ' || { echo "htseq gencode failed"; exit 1; }'."\n";
  print $sh "htseq-count --quiet --idattr=\"transcript_id\" $outdir/${f}_thout/accepted_hits_sorted.sam /nfs/bio/db/Homo_sapien/gencode_v16_with_contaminates/gencode.v16.annotation.gtf | grep ENST > $outdir/${f}_thout/gencode.transcript" . ' || { echo "htseq gencode2 failed"; exit 1; }'."\n";

  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 tophat.$f.sh");
#  system('qsub -q cluster@ccf29,cluster@ccf30,cluster@ccf31,cluster@ccf32,cluster@ccf33,cluster@ccf34,cluster@ccf36,cluster@ccf37,cluster@ccf38 -pe pe_slots 8 '."tophat.$f.sh");
}
exit;

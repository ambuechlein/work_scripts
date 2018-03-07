#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
my $dir = shift;
my $in;
# GSF1016-polyA-ishi-control-5_S13_R1_001.fastq.gz_thout/accepted_hits.sorted.bam GSF1016-polyA-ishi-control-5_S13_R1_001.fastq.gz_thout/unmapped.bam
die "can't find accepted_hits.sorted.bam for $dir\n" unless( -e "$dir/accepted_hits.sorted.bam");
die "can't find unampped.bam for $dir\n" unless( -e "$dir/unmapped.bam");
my $pre = $dir;
$pre =~ s/^GSF\d+(\-|_)//go;
$pre =~ s/_S\d+_R1_001.fastq.gz_thout//go;
my $bam = "$dir/accepted_hits.sorted.bam";
if($bam =~ /\.bam$/o){
  open($in, "samtools view -h $bam|") || die "Could not open $bam: $!\n";
}else{
  open($in, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
}
my %mapped;
my %unmapped;
while(<$in>){
  next if($_ =~ /^\@/);
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmapped{$id}++;
    next;
  }
  $mapped{$id}++;
}
close $in;
$bam = "$dir/unmapped.bam";
if($bam =~ /\.bam$/o){
  open($in, "samtools view -h $bam|") || die "Could not open $bam: $!\n";
}else{
  open($in, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
}
while(<$in>){
  next if($_ =~ /^\@/);
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmapped{$id}++;
    next;
  }
  $mapped{$id}++;
}
close $in;
my $totmapped = scalar keys %mapped;
my $totunmapped = scalar keys %unmapped;
my $totreads = $totmapped + $totunmapped;
my $single = 0;
my $multi = 0;
foreach my $id (keys %mapped){
  die "$id also flagged as unmapped\n$mapped{$id}\n$unmapped{$id}\n" if(defined $unmapped{$id});
  if($mapped{$id} > 1){
    $multi++;
  } else {
    $single++;
  }
}
my $head = "mappedSummaryHead.tsv";
if( not -e $head){
  open(my $hout, ">$head") or die "Can't open header file: $!\n";
  print $hout "Sample\tTotalReads\tTotalMapped\tPercentMapped\tTotalUnmapped\tPercentUnmapped\tSingleMapped\tPercentSingleMapped\tMultimapped\tPercentMultimapped\tUnmapped_and_Multimapped\n";
}
my $pmapped = $totmapped/$totreads*100;
my $punmapped = $totunmapped/$totreads*100;
my $psingle = $single/$totreads*100;
my $pmulti = $multi/$totreads*100;
my $un_and_multi = $totunmapped+$multi;
open(my $out, ">$pre.mappedReadSummary.tsv") or die "Can't open output $pre.mappedReadSummary.tsv: $!\n";
print $out "$pre\t$totreads\t$totmapped\t$pmapped\t$totunmapped\t$punmapped\t$single\t$psingle\t$multi\t$pmulti\t$un_and_multi\n";
close $out;
exit;

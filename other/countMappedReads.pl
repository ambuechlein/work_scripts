#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
my $bam = shift;
my $in;
#GSF886M-B3-GABIYI_S3_R1_001.sorted.bam
my $pre = basename($bam);
$pre =~ s/\.sam.gz//go;

#if($bam =~ /\.bam$/o){
  open($in, "samtools view -h  -F 256 $bam|") || die "Could not open $bam: $!\n";
#}else{
#  open($in, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
#}
my %mapped;
my %unmapped;
my $bact = 0;
while(<$in>){
  next if($_ =~ /^\@/);
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmapped{$id}++;
    next;
  }
  $mapped{$id}++;
  $bact++ unless($chr =~ /^chr/o);
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
  print $hout "Sample\tTotalReads\tTotalMapped\tPercentMapped\tTotalUnmapped\tPercentUnmapped\tSingleMapped\tPercentSingleMapped\tMultimapped\tPercentMultimapped\tUnmapped_and_Multimapped\tBacterial\tBacterialPercent\n";
}
my $pmapped = $totmapped/$totreads*100;
my $punmapped = $totunmapped/$totreads*100;
my $psingle = $single/$totreads*100;
my $pmulti = $multi/$totreads*100;
my $un_and_multi = $totunmapped+$multi;
my $pbact = $bact/$totreads*100;
##
$un_and_multi =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$single =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$multi =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$totmapped =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$totunmapped =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$totreads =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
$bact =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
##
open(my $out, ">$pre.mappedReadSummary.tsv") or die "Can't open output $pre.mappedReadSummary.tsv: $!\n";
print $out "$pre\t$totreads\t$totmapped\t$pmapped\t$totunmapped\t$punmapped\t$single\t$psingle\t$multi\t$pmulti\t$un_and_multi\t$bact\t$pbact\n";
close $out;
exit;

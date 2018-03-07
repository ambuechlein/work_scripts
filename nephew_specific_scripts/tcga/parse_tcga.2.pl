#!/usr/bin/perl
use strict;
use warnings;

# TCGA-24-2297-01A-01R-1568-13.gene.quantification.txt

my $file = shift;

my $aliq = '/nfs/labs/nephew/tcga/ov/bcr/biotab/clin/biospecimen_aliquot_ov.txt';
open(my $alin, $aliq) or die "Can't open biospecimen_aliquot_ov.txt: $!\n";
<$alin>;
my %aliq2sample;
while(<$alin>){
  my ($sid, $aid) = split(/\t/, $_);
  $aliq2sample{$aid} = $sid;
}
close $alin;
my $biosample = '/nfs/labs/nephew/tcga/ov/bcr/biotab/clin/biospecimen_sample_ov.txt';
open(my $sin, $biosample) or die "Can't open biospecimen_sample_ov.txt: $!\n";
<$sin>;
my %sample2type;
while(<$sin>){
  chomp $_;
  my @line = split(/\t/, $_);
  $sample2type{$line[0]} = $line[14];
}
close $sin;

my $clinpat = '/nfs/labs/nephew/tcga/ov/bcr/biotab/clin/clinical_patient_ov.txt';
open(my $cin, $clinpat) or die "Can't open clinical_patient_ov.txt: $!\n";
<$cin>;
my %pinfo;
while(<$cin>){
  chomp $_;
  my @line = split(/\t/, $_);
  $pinfo{$line[0]}{clinical_stage} = $line[8];
  $pinfo{$line[0]}{days_to_death} = $line[11];
}

open(my $fin, $file) or die "Can't open input $file: $!\n";
<$fin>;
print "Gene\tPatient ID\tAnalysis ID\tRPKM\tTumor Type\tDisease Stage\tDays to Death\n";
while(<$fin>){
  chomp $_;
  my @line = split(/\t/, $_);
  my $id = $line[1];
  $id =~ s/.gene.quantification.txt//go;
  $line[1] =~ /(TCGA\-\d+\-\d+)/o;
  my $pid = $1;
  my $hg19 = 0;
  if($id =~ /hg19/){
    $hg19 = 1;
    $id =~ s/.hg19//go;
  }
  die "$line[1]\t$id not found in aliquote file\n" if(not defined $aliq2sample{$id});
  if($hg19){
    print "$line[2] (hg19)\t$pid\t$line[0]\t$line[5]\t$sample2type{$aliq2sample{$id}}\t$pinfo{$pid}{clinical_stage}\t$pinfo{$pid}{days_to_death}\n";
  } else {
    print "$line[2]\t$pid\t$line[0]\t$line[5]\t$sample2type{$aliq2sample{$id}}\t$pinfo{$pid}{clinical_stage}\t$pinfo{$pid}{days_to_death}\n";
  }
}
exit;

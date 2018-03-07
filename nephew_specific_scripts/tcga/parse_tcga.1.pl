#!/usr/bin/perl
use strict;
use warnings;

# TCGA-61-2111-01A-01R-1568-13.exon.quantification.txt
# TCGA-61-2111-01A-01R-1568-13.gene.quantification.txt
# TCGA-61-2111-01A-01R-1568-13.spljxn.quantification.txt
my $id = shift;
$id =~ /(TCGA\-\d+\-\d+)/o;
my $pid = $1;

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
print "$pid\t$pinfo{$pid}{clinical_stage}\t$pinfo{$pid}{days_to_death}\t$sample2type{$aliq2sample{$id}}\n";
exit;

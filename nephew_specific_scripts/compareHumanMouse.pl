#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "human|h=s",
                         "mouse|m=s",
                         "out|o=s",
                        );
my %human;
my %mouse;
my %unmappedH;
my %unmappedM;
my $pre = basename($options{human});
$pre =~ s/.sorted\.bam//go;
warn "Parsing human\n";
open(my $hb, "samtools view -h $options{human} |") or die "Can't open $options{human}: $!\n";
while(<$hb>){
  next if($_ =~ /^\@/o);
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmappedH{$id}++;
    next;
  }
  $human{$id}++;
}
close $hb;
warn "Parsing mouse\n";
open(my $mb, "samtools view -h $options{mouse} |") or die "Can't open $options{mouse}: $!\n";
while(<$mb>){
  next if($_ =~ /^\@/o);
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmappedM{$id}++;
    next;
  }
  $mouse{$id}++;
} 
close $mb;

warn "Count reads\n";
my $hmapped = scalar keys %human;
my $mmapped = scalar keys %mouse;

my $hunmapped = scalar keys %unmappedH;
my $htotreads = $hmapped + $hunmapped;
my $munmapped = scalar keys %unmappedM;
my $mtotreads = $mmapped + $munmapped;
die "Human and Mouse have different total reads counts: $htotreads\t$mtotreads\n" unless($htotreads == $mtotreads);

my $uniqueH;
my $uniqueM;
my $onlyH;
my $onlyM;
my $both;
foreach my $id (keys %human){
  if(not defined $mouse{$id}){
    $onlyH++;
    $uniqueH++ if($human{$id} == 1);
  } else {
    $both++;
    delete $mouse{$id};
  }
}

foreach my $id (keys %mouse){
  $onlyM++;
  $uniqueM++ if($mouse{$id} == 1);
}
if(not -e "$options{out}/human_mouse_counts_header.tsv"){
  open(my $ht, ">$options{out}/human_mouse_counts_header.tsv") or die "Can't open header: $!\n";
  print $ht "FilePrefix\tTotalReads\tHumanMapped\tOnlyMappedHuman\tUniqueMappedHuman\tPercentHuman\tPercentUniqueHuman\tMouseMapped\tOnlyMouseMapped\tUniqueMappedMouse\tPercentMouse\tPercentUniqueMouse\tMappedBothHumanAndMouse\n";
  close $ht;
}
my $phuman = $hmapped/$htotreads*100;
my $pUH = $uniqueH/$htotreads*100;

my $pmouse = $mmapped/$htotreads*100;
my $pUM = $uniqueM/$htotreads*100;
open(my $ot, ">$options{out}/$pre.human_mouse_counts.tsv") or die "Can't open output $pre.human_mouse_counts.tsv: $!\n";
print $ot "$pre\t$htotreads\t$hmapped\t$onlyH\t$uniqueH\t$phuman\t$pUH\t$mmapped\t$onlyM\t$uniqueM\t$pmouse\t$pUM\t$both\n";
close $ot;
exit;

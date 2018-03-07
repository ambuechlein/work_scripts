#!/usr/bin/env perl
use strict;
use warnings;
my $bam = shift;
my $in;

die "$bam does not exist\n" unless(-e $bam);

if($bam =~ /\.bam$/o){
  open($in, "samtools view -h $bam|") || die "Could not open $bam: $!\n";
}else{
  open($in, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
}
my %best;
my %unmapped;
while(<$in>){
  if($_ =~ /^\@/){
    print $_;
    next;
  }
  my ($id, $bit, $chr) = split(/\t/, $_);
  if($chr eq '*'){
    $unmapped{$id} = $_;
    next;
  }
  $_ =~ /AS\:i\:(\-?\d+)/o;
  warn "no score\n" unless(defined $1);
  my $score = $1;
#GATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCAT      BCBFFFFFHHHHHHIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ      AS:i:-22        XN:i:0  XM:i:1  XO:i:1  XG:i:4  NM:i:5  MD:Z:3C42       YT:Z:UU
  if(not defined $best{$id} or $score > ${$best{$id}[0]}){
    ${$best{$id}[0]} = $score;
    ${$best{$id}[1]} = $_;
  }elsif($score == ${$best{$id}[0]}){
    ${$best{$id}[1]} .= $_;
  }
}
close $in;

foreach my $id (keys %best){
  die "$id also flagged as unmapped\n${$best{$id}[1]}\n$unmapped{$id}\n" if(defined $unmapped{$id});
  print ${$best{$id}[1]};
}
foreach my $id (keys %unmapped){
  print $unmapped{$id};
}

exit;

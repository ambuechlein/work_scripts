#!/usr/bin/env perl
use strict;
use warnings;
my $gtf = shift;
open(my $in, $gtf) or die "Can't open $gtf: $!\n";
while(<$in>){
  next if($_ =~ /^\#/o);
  chomp $_;
  my @d = split(/\t/, $_);
  if ($d[2] eq "gene") {
    $d[8] =~ /gene_id .(.*?).;.*gene_type .(.*?).;.*gene_name .(.*?).;/;
    my @line;
    push @line,$d[0],$d[1],$d[2],$d[3],$d[4],".",$d[6],$d[7],"ID=$1;Note=$2;Name=$3";
    my $l = join "\t",@line;
    print $l,"\n";
  } elsif ($d[2] eq "exon") {
    $d[8] =~ /gene_id .(.*?).;.*transcript_id .(.*?).;.*gene_type .(.*?).;/;
    my @line;
    push @line,$d[0],$d[1],$d[2],$d[3],$d[4],".",$d[6],$d[7],"Parent=$2";
    my $l = join "\t",@line;
    print $l,"\n";
  } elsif ($d[2] eq "CDS") {
    $d[8] =~ /gene_id .(.*?).;.*transcript_id .(.*?).;.*gene_type .(.*?).;/;
    my @line;
    push @line,$d[0],$d[1],$d[2],$d[3],$d[4],".",$d[6],$d[7],"Parent=$2";
    my $l = join "\t",@line;
    print $l,"\n";
  } elsif ($d[2] eq "transcript") {
    $d[8] =~ /gene_id .(.*?).;.*transcript_id .(.*?).;.*gene_type .(.*?).;/;
    my @line;
    push @line,$d[0],$d[1],"mRNA",$d[3],$d[4],".",$d[6],$d[7],"ID=$2;Parent=$1;Name=$2.$3";
    my $l = join "\t",@line;
    print $l,"\n";
  }
}

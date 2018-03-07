#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
my $bam = shift;
my $in;
#GSF886M-B3-GABIYI_S3_R1_001.sorted.bam
my $pre = basename($bam);
$pre =~ s/\.sorted\.bam//go;

if($bam =~ /\.bam$/o){
  open($in, "samtools view -h $bam|") || die "Could not open $bam: $!\n";
}else{
  open($in, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
}
my %mapped;
my %read;
while(<$in>){
  if($_ =~ /^\@/){
    print $_;
    next;
  }
  my ($id, $bit, $chr) = split(/\t/, $_);
  next if($chr eq '*');
  $mapped{$id}++;
#  $read{$id} = $_;
  push(@{$read{$id}},$_);
}
close $in;

foreach my $id (keys %mapped){
  if($mapped{$id} <= 1){
    print ${$read{$id}}[0];
  }else{
    my $size = scalar @{$read{$id}};
    my $idx = int(rand($size));
    print ${$read{$id}}[$idx];
  }
}
exit;

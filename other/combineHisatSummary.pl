#!/usr/bin/env perl
use strict;
use warnings;

my $dir = shift;
my $so = shift;
my @order;
open(my $sin, $so) or die "Can't open $so: $!\n";
while(<$sin>){
  chomp $_;
  my ($x,$y) = split(/\t/,$_);
  $x .= '.paired.fastq.gz_hisat';
  push(@order, $x);
}
opendir(my $din, $dir) or die "Can't open directory $dir: $!\n";
my %map;
my @fo;
while(my $pre = readdir($din)){

  next unless($pre =~ /.paired.fastq.gz_hisat$/go);
#  next unless($pre =~ /001.fastq.gz_hisat$/go);
  my $f = $pre;
 $f =~ s/_R1_001.paired.fastq.gz_hisat$/\.summary\.tsv/go;
#  $f =~ s/_R1_001.fastq.gz_hisat$/\.summary\.tsv/go;
  die "no summary $pre/$f\n" unless(-e "$pre/$f");
  open(my $in, "$pre/$f") or die "Can't open $pre/$f: $!\n";
  push(@fo,$pre);
  <$in>;
  while(<$in>){
    chomp $_;
    $_ =~ s/^\s+//go;
    my ($t,$c) = split(/\: /,$_);
    $map{$t}{$pre} = $c;
  }
  close $in;
}
closedir $din;

foreach my $pre (@order){
  print "\t$pre";
}
print "\n";
my @mapOrder = ('Total pairs',
                'Aligned concordantly 1 time',
                'Aligned concordantly >1 times',
                'Aligned concordantly or discordantly 0 time',
                'Aligned discordantly 1 time',
                'Total unpaired reads',
                'Aligned 0 time',
                'Aligned 1 time',
                'Aligned >1 times',
                'Overall alignment rate');
#  foreach my $t (keys %map){
foreach my $t (@mapOrder){
  print $t;
  foreach my $pre (@order){
    print "\t$map{$t}{$pre}";
  }
  print "\n";
}
exit;

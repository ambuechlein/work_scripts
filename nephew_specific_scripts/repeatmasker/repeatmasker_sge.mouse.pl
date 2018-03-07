#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $results = GetOptions(\%options,
                         "list|l=s",
                         "outdir|o=s",
                        );
my $list = $options{list};
open(my $in, $list) or die "Can't open list file: $!\n";

while(<$in>){
  chomp $_;
  my $command = "/home/abuechle/bin/RepeatMasker/RepeatMasker -e crossmatch -pa 4 -species mouse -no_is -frag 300000 -dir $options{outdir} -small -xm $_";
#  $_ =~ /\/nfs\/labs\/nephew\/6-13-13\/repeatMasker\/reads\/(\S+).nonredundant.fasta(\d+).fsa/o;
#  my $base = 'repeats_' . $1.'_'.$2;
  my $base = basename($_);
  $base = 'repeats_' . $base;
  open(my $sh, ">$base.sh") or die "Can't open shell script $base.sh: $!\n";
  print $sh '#!/bin/bash', "\n\n";
  print $sh "$command". ' || { echo "RepeatMasker failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;
  system("chmod 755 $base.sh");
  system("qsub -q cluster -p -1000 -pe pe_slots 4 -wd $cwd $base.sh");
}
close $in;
exit;

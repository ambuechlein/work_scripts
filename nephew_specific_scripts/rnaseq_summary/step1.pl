#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $results = GetOptions(\%options,
                         "dir|d=s",
                         "gtf|g=s",
                         "tophat_dir|t=s",
                         "chr_list|c=s",
                        );
my $workingDir = $options{dir};
my $gtfFile = $options{gtf};
system("mkdir -p $workingDir/result");
system("find $options{tophat_dir}/* -follow -name accepted_hits_sorted.bam > $workingDir/bam.list");

open(my $cin, $options{chr_list}) or die "Can't open chromosome list $options{chr_list}: $!\n";
my @file;
my $c;
my $fh;
while (<$cin>) { 
  chomp $_; 
  my $chr = $_;
  $c++;
  $c = $c%26;
  open($fh, ">map.$chr.sh") or die "Can't open map.$c.sh: $!\n";
  print $fh '#!/bin/bash', "\n";
  print $fh 'echo ">>>>> step1 tophat2gff "', "\n";
  print $fh 'echo ">>>>> startdate "`date`', "\n";
  print $fh 'echo ">>>>> hostname "`hostname`', "\n";
  print $fh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $fh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $fh 'echo "Starting mapStrandedTopHat2GFF5.pl"', "\n";
  my $let = chr(65+$c); 
  my $cmd = "perl /nfs/labs/nephew/scripts/rnaseq_summary/mapStrandedTopHat2GFF5.new.pl $gtfFile $chr $workingDir/bam.list $workingDir/result/ $let$let" . ' || { echo "Step1 tophat2gff mapping  failed"; exit 1; }'."\n";
  print $fh $cmd,"\n";
  print $fh 'echo "Finished"', "\n";
  system("chmod 755 map.$chr.sh");
  if($chr eq 'chr1'){
    system('qsub -q bigmem@exon,bigmem@intron -pe pe_slots 8 '."-wd $cwd map.$chr.sh");
  } elsif($chr =~ /^chr/o){
    system('qsub -q bigmem@exon,bigmem@intron -pe pe_slots 4 '."-wd $cwd map.$chr.sh");
  } else {
    system("qsub -q cluster,bigmem -pe pe_slots 8 -wd $cwd map.$chr.sh");
  }
}

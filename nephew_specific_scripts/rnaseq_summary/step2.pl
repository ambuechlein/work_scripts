#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $cwd = cwd();

my %options;
my $results = GetOptions(\%options,
                         "dir|d=s",
                         "samplelist|s=s",
                        );
my $dir = $options{dir};
open(my $sin, $options{samplelist}) or die "Can't open sample list $options{samplelist}: $!\n";
my @file;
my $c;
my $fh;
while (<$sin>) { 
# for i in `cat sample.list`; do echo $i; echo -e '#!/bin/bash\n\n' > sort.$i.sh; echo "sort -u -S 2% -T . /nfs/labs/nephew/burow/LKB1/summary/result/$i.*.counts | gzip -c > /nfs/labs/nephew/burow/LKB1/summary/result/$i.counts.sorted.gz" >> sort.$i.sh; done;
  chomp $_; 
  my $i = $_;
  open($fh, ">sort.$i.sh") or die "Can't open map.$c.sh: $!\n";
  print $fh '#!/bin/bash', "\n";
  print $fh 'echo ">>>>> step2 sort "', "\n";
  print $fh 'echo ">>>>> startdate "`date`', "\n";
  print $fh 'echo ">>>>> hostname "`hostname`', "\n";
  print $fh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $fh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $fh 'echo "Starting Sort"', "\n";
  my $cmd = "sort --parallel 8 -u -S 10% -T . $dir/$i.*counts | sort --parallel 8 -k 1,1 -S 10% -T . | gzip -c > $dir/$i.counts.sorted.gz" . ' || { echo "Step2 sorting counts  failed"; exit 1; }'."\n";

  print $fh $cmd,"\n";
  print $fh 'echo "Finished"', "\n";
  system("chmod 755 sort.$i.sh");
#  system('qsub -q cluster,bigmem@antigen,bigmem@antibody,bigmem@trna -pe pe_slots 8 '."sort.$i.sh");
  system('qsub -q cluster -pe pe_slots 8 '." -wd $cwd sort.$i.sh");
}
close $sin;
exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Data::Dumper;

BEGIN {
  my $cwd = cwd();
  $ENV{'SGE_O_HOME'} = $cwd;
};
my %options;
my $result = GetOptions(\%options,
                        "control|c=s",
                        "table|t=s",
                        "pheno|p=s",
                       );

my $cwd = cwd();
my @de = ("newDESeq2.R", "newDESeq.R", "newEdgeR.R");
my $jobStates;
foreach my $sh (@de){
  my $cmd = &createSGEScript($sh, $options{control}, $options{table}, $options{pheno});
  my $sto = `qsub -q bigmem -pe pe_slots 16 -wd $cwd $cwd/$cmd`;
  chomp $sto;
  $sto =~ /Your job (\d+) .* has been submitted/io;
  my $id = $1;
##  $jobStates->{'job'}->{$id} = 'submitted'
  $jobStates->{$id}->{status} = 'submitted'
}
my $poll = 0;
while($poll < 3){
  &pollJob($jobStates);
#  print Dumper($jobStates), "\n";
  $poll = 0;
  foreach my $pid (keys %{$jobStates}){
    $poll++ if($jobStates->{$pid}->{status} eq 'finished');
    print "$pid\t$jobStates->{$pid}->{status}\n" if($jobStates->{$pid}->{status} ne 'finished');
  }
  print "\n";
  next if($poll > 2);
  sleep(10);
}

$poll = 0;
foreach my $pid (keys %{$jobStates}){
  my $file = $jobStates->{$pid}->{logs} .".o$pid";
  my $line = `tail -n 1 $file`;
  chomp $line;
  if($line eq 'Finished'){
    $poll++;
  }else{
    warn "Please investigate $file\nJob $pid possibly failed\n";
  }
}
if($poll == 3){
  print "All Jobs Done\n";
} else {
  print "Possible job failures\n";
  exit(1);
}
exit;

sub pollJob{
  my ($jobStates) = @_;
#  foreach my $pid (keys %{$jobStates->{'job'}}){
  foreach my $pid (keys %{$jobStates}){
    my @status = `qstat -j $pid 2>&1`;
    while (@status) {
      my $line=shift @status;
      chomp($line);
      next if($line =~ /^\=+$/o); 
      if($line =~ /Your job 4843295 \(\.+\) has been submitted/){
#        $jobStates->{'job'}->{$pid}='pending'
        $jobStates->{$pid}->{status}='pending'
      }elsif($line =~ /Following jobs do not exist:/){
        while(@status){
          my $line2 = shift @status;
          chomp $line2;
          die "We are polling the wrong ID somehow: $pid $line2" unless($pid == $line2);
#          $jobStates->{'job'}->{$pid}='finished'
          $jobStates->{$pid}->{status}='finished'
        }
      }elsif ($line =~ /^job_number\:\s+(\d+)/) {
        die "We are polling the wrong ID somehow: $pid $1" unless($pid == $1);
#        $jobStates->{'job'}->{$pid}='running';
        $jobStates->{$pid}->{status}='running';
        my $logs = '';
        while(@status){
          my $line2 = shift @status;
          chomp $line2;
          if($line2 =~ /cwd\:\s+(\S+)/o){
             $logs = $1;
          }elsif($line2 =~ /job_name\:\s+(\S+)/){
             $logs .= "/$1";
          }
#       submission_time:            Fri Dec 12 10:03:54 2014
#       owner:                      abuechle
#       uid:                        2907
#       group:                      cgb
#       gid:                        1000
#       sge_o_home:                 /home/abuechle
#       sge_o_log_name:             abuechle
#       sge_o_path:                 /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/home/abuechle/bin/cutadapt-1.4.1/bin:/nfs/bio/sw/encap/bowtie-0.12.8/bin:/nfs/bio/sw/bin:/home/abuechle/bin:/home/abuechle/bin/miRDP1.3:/home/abuechle/bin/ViennaRNA/bin:/home/abuechle/.gsh:/nfs/projects/isga/vtest/workflow-3.2.1-cgbr0/apache-ant-1.8.1/bin:/usr/local/bin:/cluster/solexa/bin:/home/abuechle/.local/bin:/home/abuechle/bin/moabs-v1.2.9.src.x86_64_Linux.data/bin
#       sge_o_shell:                /bin/tcsh
#       sge_o_workdir:              /nfs/projects/bioinfo/abuechle
#       sge_o_host:                 antigen
#       account:                    sge
#       mail_list:                  abuechle@antigen
#       notify:                     FALSE
#       job_name:                   date.sh
#       priority:                   -100
#       jobshare:                   0
#       hard_queue_list:            bigmem
#       env_list:                   
#       project:                    global
#       scheduling info:
        }
        $jobStates->{$pid}->{logs}=$logs;
      }
    }
  }
}

sub createSGEScript {
  my ($cmd, $ctl, $tbl, $pheno) = @_;
#  my $cwd = cwd()
  open(my $sh, ">$cmd.SGE.$$.sh") or die "can't open shell script $cmd.SGE.$$.sh: $!\n\n";

  print $sh '#!/bin/bash';
  print $sh "\n";
  print $sh 'echo ">>>>> '.$cmd.' "', "\n";
  print $sh 'echo ">>>>> startdate "`date`', "\n";
  print $sh 'echo ">>>>> hostname "`hostname`', "\n";
  print $sh 'echo ">>>>> job_name "$JOB_NAME', "\n";
  print $sh 'echo ">>>>> job_id "$JOB_ID', "\n";

  print $sh 'echo "Starting '.$cmd.'"', "\n";
  print $sh "Rscript /home/abuechle/bin/$cmd $ctl $tbl $pheno". ' || { echo "'.$cmd.' failed"; exit 1; }'."\n";
  print $sh 'echo "Finished"', "\n";
  close $sh;

  system("chmod 755 $cmd.SGE.$$.sh");
  return "$cmd.SGE.$$.sh";
}

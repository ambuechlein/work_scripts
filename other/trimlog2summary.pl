#!/usr/bin/env perl
use strict;
use warnings;
use Number::Format 'format_number';

my %map;
print "Sample\tFile\tInputReadPairs\tBothSurviving\tBothSurviving%\tForwardOnlySurviving\tForwardOnlySurviving%\tReverseOnlySurviving\tReverseOnlySurviving%\tDropped\tDropped%\n";

foreach my $f (@ARGV){
  open(my $fin, $f) or die "Can't open $f: $!\n";
  my $pre = $f;
# trim.GSF1429-BC-2_S4_R1_001.fastq.gz.sh.log
  $pre =~ s/trim\.//go;
  $pre =~ s/_R1_001.fastq.gz.sh.*//go;
  my $x = $pre;
  $x =~ s/GSF\d+\-//g;
  $x =~ s/_S\d+//g;
  while(<$fin>){
    chomp $_;
#   trim.GSF1429-BC-2_S4_R1_001.fastq.gz.sh.log:Input Read Pairs: 30372797 Both Surviving: 29288062 (96.43%) Forward Only Surviving: 979931 (3.23%) Reverse Only Surviving: 87922 (0.29%) Dropped: 16882 (0.06%)
    next unless($_ =~ /^Input Read/g);
    $_ =~ /Input Read Pairs: (\d+) Both Surviving: (\d+) \((\S+)\) Forward Only Surviving: (\d+) \((\S+)\) Reverse Only Surviving: (\d+) \((\S+)\) Dropped: (\d+) \((\S+)\)/;

    # $map{$pre}{InputReadPairs} = $1;
    # $map{$pre}{BothSurviving} = $2;
    # $map{$pre}{'BothSurviving%'} = $3;
    # $map{$pre}{ForwardOnlySurviving} = $4;
    # $map{$pre}{'ForwardOnlySurviving%'} = $5;
    # $map{$pre}{ReverseOnlySurviving} = $6;
    # $map{$pre}{'ReverseOnlySurviving%'} = $7;
    # $map{$pre}{Dropped} = $8;
    # $map{$pre}{'Dropped%'} = $9;
#    print "$x\t$pre\t$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\n";
    print "$x\t$pre"; #$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\n";
    print "\t",format_number($1)
    print "\t",format_number($2)
    print "\t",format_number($3)
    print "\t",format_number($4)
    print "\t",format_number($5)
    print "\t",format_number($6)
    print "\t",format_number($7)
    print "\t",format_number($8)
    print "\t",format_number($9),"\n";

  }
  close $fin;
}
exit;

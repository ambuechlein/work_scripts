#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
#####################################################################################################################################################
#  C Braun      CGB m/d/y
#  Name         domethylTable.pl
#  Purpose      Take in methylation call file, use GFF to define genes and describe methylation in body, up and down
#  Input        ITERATIVELY - ONE CHR AT A TIME!!!!
#  Output       
#  Sample       perl doMethylTable.pl -c chr1 -i /nfs/labs/nephew/11-25-13-Run/methylation/bismark/131125_RKN01B3_1_TATAAT_L002_R1/CpG_context_131125_RKN01B3_1_TATAAT_L002_R1.fastq_bismark_bt2.txt -o output -g /nfs/bio/db/Homo_sapien/gencode_v19/jbrowse.gencode.v19.annotation.gff -t CpG -h yes
#####################################################################################################################################################

my %options;
my $result = GetOptions(\%options,
                        "chr|c=s",
                        "input|i=s",
                        "output|o=s",
                        "gtf|g=s",
                        "type|t=s",
                        "header|h=s",
                       );
my $chr = $options{chr};
my $type = $options{type};
my $reference = $options{gtf};
my $input = $options{input};
my $output = $options{output}.".$chr.$type.tsv";

my ($GeneID, $GeneStart, $GeneEnd, $GeneStrand, $GeneChr, @GeneArray);
my $count = 0;
open (my $out, ">$output") or die "Can't open $output for writing: $!\n";
print $out "GeneID\tGeneChr\tGeneStart\tGeneEnd\tGeneStrand\tBodyC\tBodyMet\tBodyProp\tUpC\tUpMet\tUpProp\tDownC\tDownMet\tDownProp\n" if ($options{header} eq "yes");
print "Loading BRAT files\n";
my $MethylHash = {};

open (my $in, "<$input") or die ("Could not open input file $input: $!\n");
<$in>;
while (my $line = (<$in>)){
   chomp($line);
   my @temp = split /\t/, $line;
#  HWI-ST1052:152:D2GYAACXX:2:1101:1193:2245_1:N:0:TATAAT  +       chr9    32689014        Z
   next unless($temp[2] eq $chr);
     ${$MethylHash->{$temp[1]}->{$temp[3]}}[1]++;
     ${$MethylHash->{$temp[1]}->{$temp[3]}}[0] = 0 if(not defined ${$MethylHash->{$temp[1]}->{$temp[3]}}[0]);
     if ( $temp[4] eq "\U$temp[4]" ){ #is Methylated
       ${$MethylHash->{$temp[1]}->{$temp[3]}}[0]++;
     }
}
close ($in);
print "Input loaded\n";

################ Open the reference file, read just the relevant chromosome and pull data from hash

open (my $ein, "<$reference") or die ("Could not open reference file, sorry");
   while (my $line = (<$ein>)){
      chomp($line);
      my @inside = split /\t/, $line;
      next if ($inside[0] ne $chr);  ### short circuit if not on correct chromosome

      if ($inside[2] eq "gene"){
         $inside[8] =~ /ID\=(\w+\d+\.\d+)/;  ### grab the geneID
         $GeneID = $1;
         $GeneStart = $inside[3];
         $GeneEnd = $inside[4];
         $GeneStrand = $inside[6];
         $GeneChr = $inside[0];
         countMethyl($MethylHash, $GeneID, $GeneChr, $GeneStart, $GeneEnd, $GeneStrand);
      }
   }
close ($out);
close ($ein);
exit;

sub countMethyl {
  my ($MethylHash, $GeneID, $GeneChr, $GeneStart, $GeneEnd, $GeneStrand) = @_;
  my ($BCcount, $BMcount, $BNcount, $BProp, $UCcount, $UMcount, $UNcount, $UProp, $DCcount, $DMcount, $DNcount, $DProp) = (0) x 12;
  my $place = ($GeneStart - 500);
  my $DownEnd = ($GeneEnd + 500);
  while ($place < $GeneStart){
     if ($MethylHash->{$GeneStrand}->{$place}){
        $UCcount++;
        my $percent = ${$MethylHash->{$GeneStrand}->{$place}}[0] / ${$MethylHash->{$GeneStrand}->{$place}}[1];
        if ($percent >= 0.5){
           $UMcount++;
        }else{
           $UNcount++;
        }
     }
     $place++;
  }
  if ($UCcount == 0){
     $UProp = 0;
  }else{
     $UProp = ($UMcount / $UCcount);
  }

  while ($place <= $GeneEnd){
     if ($MethylHash->{$GeneStrand}->{$place}){
        $BCcount++;
        my $percent = ${$MethylHash->{$GeneStrand}->{$place}}[0] / ${$MethylHash->{$GeneStrand}->{$place}}[1];
        if ($percent >= 0.5){
           $BMcount++;
        }else{
           $BNcount++;
        }
     }
     $place++;
  }
  if ($BCcount == 0){
     $BProp = 0;
  }else{
     $BProp = ($BMcount / $BCcount);
  }

  while ($place <= $DownEnd){
     if ($MethylHash->{$GeneStrand}->{$place}){
        $DCcount++;
        my $percent = ${$MethylHash->{$GeneStrand}->{$place}}[0] / ${$MethylHash->{$GeneStrand}->{$place}}[1];
        if ($percent >= 0.5){
           $DMcount++;
        }else{
           $DNcount++;
        }
     }
     $place++;
  }
  if ($DCcount == 0){
     $DProp = 0;
  }else{
     $DProp = ($DMcount / $DCcount);
  }
  print $out "$GeneID\t$GeneChr\t$GeneStart\t$GeneEnd\t$GeneStrand\t$BCcount\t$BMcount\t$BProp\t$UCcount\t$UMcount\t$UProp\t$DCcount\t$DMcount\t$DProp\n";
}

#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

#####################################################################################################################################################
#  C Braun      CGB m/d/y
#  Name         domethylTable.pl
#  Purpose      Take in methylation call file, use GFF to define genes and describe methylation in body, up and down
#  Input        ITERATIVELY - ONE CHR AT A TIME!!!!
#  Output       
#  Sample       ./doMethylTable.pl -s chr1 -f ~/cb_directory/bisulfite_processing/BRAT/nephew/BT04/results_dir/prefix_new_forw.txt -r ~/cb_directory/bisulfite_processing/BRAT/nephew/BT04/results_dir/prefix_new_rev.txt -o output -g /nfs/labs/nephew/drusch/nextAnalysis/result/xxx.gtf -t CHH 
#####################################################################################################################################################

my %options;
my $result = GetOptions(\%options,
                        "chr|c=s",
                        "forward|f=s",
                        "reverse|r=s",
                        "output|o=s",
                        "gtf|g=s",
                        "type|t=s",
                        "header|h=s",
                       );
my $chr = $options{chr};
my $header = $options{header};
my $type = $options{type};
my $reference = $options{gtf};
my $forward = $options{forward};
my $reverse = $options{reverse};
my $output = $options{output}.".$chr.$type.tsv";

my ($GeneID, $GeneStart, $GeneEnd, $GeneStrand, $GeneChr, @GeneArray);
my $count = 0;
open (my $out, ">$output") or die "Can't open $output for writing: $!\n";

if ($header eq "yes"){
   print $out "GeneID\tGeneChr\tGeneStart\tGeneEnd\tGeneStrand\tBodyC\tBodyMet\tBodyProp\tUpC\tUpMet\tUpProp\tDownC\tDownMet\tDownProp\n";
}
my $MethylHash = {};
print "Loading BRAT files\n";
parseBRAT($MethylHash, $forward, $chr, $type);
print "Forward strand loaded\n";
parseBRAT($MethylHash, $reverse, $chr, $type);
print "Reverse strand loaded\n";

################ 3) Open the reference file, read just the relevant chromosome and pull data from hash

open (my $ein, "<$reference") or die ("Could not open reference file, sorry");
while (my $line = (<$ein>)){
   chomp($line);
   my @inside = split /\t/, $line;
   next if ($inside[0] ne $chr);  ### short circuit if not on correct chromosome
   next unless ($inside[2] eq "gene");
   $inside[8] =~ /ID\=(\w+\d+\.\d+)/;  ### grab the geneID
   $GeneID = $1;
   $GeneStart = $inside[3];
   $GeneEnd = $inside[4];
   $GeneStrand = $inside[6];
   $GeneChr = $inside[0];
   my ($BCcount, $BMcount, $BNcount, $BProp, $UCcount, $UMcount, $UNcount, $UProp, $DCcount, $DMcount, $DNcount, $DProp) = (0) x 12;
   my $place = ($GeneStart - 500);
   my $DownEnd = ($GeneEnd + 500);   
   while ($place < $GeneStart){
      if ($MethylHash->{$GeneStrand}->{$place}){
         $UCcount++;
         my $now = $MethylHash->{$GeneStrand}->{$place};
         $now =~ /(\d+):(\S+)/; #die "NOW 196: $now\n";
         if ($2 >= 0.5){
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
         my $now = $MethylHash->{$GeneStrand}->{$place};
         $now =~ /(\d+):(\S+)/;
         if ($2 >= 0.5){
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
         my $now = $MethylHash->{$GeneStrand}->{$place};
         $now =~ /(\d+):(\S+)/;
         if ($2 >= 0.5){
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
close ($out);
close ($ein);
exit;

sub parseBRAT{
  my ($MethylHash, $forward, $chr, $type) = @_;
  my $streak = "no";
  open (my $fin, "<$forward") or die ("Could not open forward file, sorry");

  while (my $line = (<$fin>)){
     chomp($line);
     my @temp = split /\t/, $line;
# chr1    10003   10003   CHH:0   0       + 
     if (($streak eq "yes") and ($temp[0] ne $chr)){
        print "\nEnd of Streak\n";
        last;
     }
  
     if ($temp[0] eq $chr){
        $streak = "yes";
        $temp[3] =~ /(\w+):(\d+)/;
        my $kind = $1;
        my $Cbase = $2;
        my $Ctotal = $temp[4];
  
        if ($type eq $kind){
           if (($Cbase != 0) or ($Ctotal != 0)){
              my $str = $Cbase.":".$Ctotal;
              $MethylHash->{$temp[5]}->{$temp[1]} = $str;
           }
        }
     }
  }
  close ($fin);
#  return;
}

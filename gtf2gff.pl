#!/usr/bin/env perl
use strict;
use warnings;
my $inG = shift;
my $outG = shift;
open(my $in, $inG) or die "Can't open file $inG: $!\n";
while(<$in>){
  chomp $_;
# chr1	hg19_rmsk	exon	16777161	16777470	2147.000000	+	.	gene_id "AluSp"; transcript_id "AluSp"; 
# chr1	hg19_rmsk	exon	25165801	25166089	2626.000000	-	.	gene_id "AluY"; transcript_id "AluY"; 
# chr1	hg19_rmsk	exon	33553607	33554646	626.000000	+	.	gene_id "L2b"; transcript_id "L2b";
# chr1	HAVANA	gene	11869	14412	.	+	.	ID=ENSG00000223972.4;Note=pseudogene;Name=DDX11L1

  my @d = split(/\t/, $_);
  if ($d[2] eq "exon") { 
    $d[8] =~ /gene_id .(.*?).;.*transcript_id .(.*?).;/; 
    my @l = ();
    push @l,$d[0],$d[1],$d[2],$d[3],$d[4],".",$d[6],$d[7],"ID=$1;NAME=$2;"; 
    my $pl = join "\t",@l; 
    print $pl,"\n"; 
  }
}
close $in;
exit;

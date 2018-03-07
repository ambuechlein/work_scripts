#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my %options;
my $results = GetOptions ( \%options,
                           'input|i=s'
                         );
my $infile = $options{input};
my $infileformat = 'fasta';
# >comp0_c0_seq1 len=204 path=[1:0-203]

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my %max;
my %lengths;
# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
#   print $inseq->id, "\t", $inseq->length(), "\n";
    my $id = $inseq->id;
    $id =~ /(comp\d+_c\d+)_seq/o;
    my $isogroup = $1;
    if( (not defined $lengths{$isogroup}) or ($inseq->length() > $lengths{$isogroup}) ){
      $lengths{$isogroup} = $inseq->length();
      $max{$isogroup} = $inseq;
    }
}

my $outfile = $infile . '.new';
my $outfileformat = 'fasta';
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                              '-format' => $outfileformat);
use Data::Dumper;
foreach my $inseq (sort keys %max){
#   warn Dumper($max{$inseq});
   $seq_out->write_seq($max{$inseq});
}
exit;

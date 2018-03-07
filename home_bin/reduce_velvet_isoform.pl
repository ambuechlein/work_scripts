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
# Locus_4041_Transcript_9/10_Confidence_0.044_Length_1734
# Locus_4041_Transcript_10/10_Confidence_0.400_Length_3124
# Locus_4042_Transcript_1/4_Confidence_0.667_Length_2427
# Locus_4042_Transcript_2/4_Confidence_0.444_Length_1393


# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my %max;
my %lengths;
# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
#   print $inseq->id, "\t", $inseq->length(), "\n";
    my $id = $inseq->id;
    $id =~ /(Locus_\d+)_Transcript_\d+\/\d+/o;
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

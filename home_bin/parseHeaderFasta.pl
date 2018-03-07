#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
# get command-line arguments, or die with a usage statement
use Getopt::Long;
my %options;
my $results = GetOptions ( \%options,
                           'input|i=s'
                         );
die "Please provide and input file\nUsage: perl refromatFasta.pl -i <fasta_file>\n\n" unless $options{input};
my $infile = $options{input};
my $infileformat = 'fasta';
my $outfile = $infile . '.new';
my $outfileformat = 'fasta';

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                              '-format' => $outfileformat);

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
   $seq_out->write_seq($inseq);
}
exit;

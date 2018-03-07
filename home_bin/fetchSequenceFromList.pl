#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
# get command-line arguments, or die with a usage statement
use Getopt::Long;
my %options;
my $results = GetOptions ( \%options,
                           'id|i=s',
                           'fasta|f=s'
                         );
open(my $in, $options{id}) or die "can't open id file: $!\n";
my %ids;
while(<$in>){
  chomp $_;
  $ids{$_} = 1;
}
my $infile = $options{fasta};
my $infileformat = 'fasta';
my $outfile = 'unannotated.fasta';
my $outfileformat = 'fasta';

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                              '-format' => $outfileformat);

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
   my $x = $inseq->primary_id();
   $seq_out->write_seq($inseq) if(defined $ids{$x} );
}
exit;

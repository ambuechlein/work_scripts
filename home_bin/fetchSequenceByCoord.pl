#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
# get command-line arguments, or die with a usage statement
use Getopt::Long;
use Data::Dumper;
my %options;
my $results = GetOptions ( \%options,
                           'id|i=s',
                           'fasta|f=s',
                           'start|s=s',
                           'end|e=s',
                           'outname|o=s',
                           'strand|d=s',
                         );
die "Please provide and input file\nUsage: perl refromatFasta.pl -i <fasta_file>\n\n" unless $options{fasta};
my $infile = $options{fasta};
my $infileformat = 'fasta';
my $outfile = $options{outname} . '.fasta';
my $outfileformat = 'fasta';

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                              '-format' => $outfileformat);

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
#   warn $inseq->primary_id(), "\n";
   if($inseq->primary_id() eq $options{id}){
     warn $inseq->primary_id(), "\n";
     my $substring = $inseq->subseq($options{start},$options{end});

     my $seqobj = Bio::PrimarySeq->new (
         -seq              => $substring,
         -id               => $options{outname},
         -alphabet         => 'dna',
         -is_circular      => 0,
     );
     if($options{strand} eq '+'){
       $seq_out->write_seq($seqobj);
     } else {
       $seq_out->write_seq($seqobj->revcom());
     }
   }
}
exit;

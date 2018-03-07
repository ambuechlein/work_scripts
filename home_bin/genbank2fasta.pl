#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
my $infile  = shift;
my $outfile   = shift;
 
my $seqin = Bio::SeqIO->new(
                            -file     => $infile,
                            -format => "genbank",
                            );
my $seqout = Bio::SeqIO->new(
                             -file   => ">$outfile",
                             -format => "fasta",
                             );
 
my $seq_object = $seqin->next_seq;
 
for my $feat_object ($seq_object->get_SeqFeatures) {
   if ($feat_object->primary_tag eq "CDS") {
        if($feat_object->has_tag('translation')){
            my ($locus_tag, $gene, $product) = '';
            my $start = $feat_object->location->start;
            my $end = $feat_object->location->end;
            my $strand = $feat_object->location->strand == 1 ? "$start:$end Forward" : "$start:$end Reverse";
            for my $val ($feat_object->get_tag_values('locus_tag')){ $locus_tag = $val };
            for my $val ($feat_object->get_tag_values('gene')){ $gene = $val };
            for my $val ($feat_object->get_tag_values('product')){ $product = $val };
            $gene = $locus_tag if $gene eq '';
            my $seq_obj = Bio::Seq->new(-seq => $feat_object->get_tag_values('translation'),
                                        -display_id => $locus_tag.' | '.$gene.' | '.$product.' | '. $strand );
            $seqout->write_seq($seq_obj);
        }
   }
}

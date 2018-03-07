#!/usr/bin/env perl
use strict;
use warnings;
use lib "/data/web/isga.cgb/docs/lib/perl5/";
use ISGA;
use Data::Dumper;

  my $nuc_ref = [];
  my $prot_ref = [];
  foreach ( @{ISGA::ReferenceDB->query( )} ){
    if($_->getTemplate->getFormat eq 'BLAST Nucleotide Database' and $_->getRelease->getStatus->getName eq 'Published'){
        push(@$nuc_ref, $_);
    } elsif($_->getTemplate->getFormat eq 'BLAST Amino Acid Database' and $_->getRelease->getStatus->getName eq 'Published') {
        push(@$prot_ref, $_);
    }
  }


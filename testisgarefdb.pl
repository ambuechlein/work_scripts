#!/usr/bin/env perl
use strict;
use warnings;
use lib "/data/web/isga.cgb/docs/lib/perl5/";
use ISGA;
use Data::Dumper;

my $ref = ISGA::Reference->new(Id => '12');
my $template = ISGA::ReferenceTemplate->new( Id => 35 );

if ( my $db = $ref->getReferencePath($template) ) {
  print Dumper($db);
} else {
  print "wahwah\n";
}

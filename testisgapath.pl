#!/usr/bin/env perl
use strict;
use warnings;
use lib "/data/web/isga.cgb/docs/lib/perl5/";
use ISGA;
my $self = ISGA::Pipeline->new(Id => '3');
print $self->getStatus->getName, "\n";

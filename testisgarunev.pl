#!/usr/bin/env perl
use strict;
use warnings;
use lib "/data/web/isga.cgb/docs/lib/perl5/";
use ISGA;
my $self = ISGA::RunEvidenceDownload->new(Id => '95');
my $path = $self->getFilePath;
print "$path\n";

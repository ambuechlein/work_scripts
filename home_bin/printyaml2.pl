#!/usr/bin/env perl
use strict;
use warnings;
use YAML;
my $analysis = shift;

use Data::Dumper;
my $data = YAML::LoadFile($analysis);
print $data->{Analysis}->{name1}->{VALUE}, "\n";
print Dumper($data);

#!/usr/bin/env perl
use strict;
use warnings;
use YAML;
local $YAML::UseHeader = 0;
my $data = {};
$$data{RunID} = 1;
#my $hash = {"NAME" => "name1", "VALUE" => "8", "DESCRIPTION" => "some long words"};
#push(@{$$data{Analysis}}, $hash);
#
#my $hash2 = {"NAME" => "name2", "VALUE" => "10", "DESCRIPTION" => "some longer words"};
#push(@{$$data{Analysis}}, $hash2);

my $hash = {};
$hash->{name1} = {"VALUE" => "8", "DESCRIPTION" => "some long words"};
#push(@{$$data{Analysis}{name1}}, $hash);
#$data->{Analysis}->{name1} = $hash;

my $hash2 = {};
$hash2->{name2} = {"VALUE" => "10", "DESCRIPTION" => "some longer words"};
#push(@{$$data{Analysis}{name2}}, $hash2);
#$data->{Analysis}->{name2} = $hash2;

foreach my $k (keys %$hash) { 
  $data->{Analysis}->{$k} = $hash->{$k};
}

foreach my $k (keys %$hash2) { 
  $data->{Analysis}->{$k} = $hash2->{$k};
}

open(my $out, ">test.yaml") or die "Can't open yaml output\n";
print $out '--- !perl/ISGA::RunAnalysis'."\n";
print $out YAML::Dump($data);

use Data::Dumper;
print Dumper($data);
close $out;

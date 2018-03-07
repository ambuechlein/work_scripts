#!/usr/bin/env perl

use File::Find;
my @files;

find(
{
wanted => sub {
push @files, $File::Find::fullname
if -f $File::Find::fullname && /\.pm$/
},
follow => 1,
follow_skip => 2,
},
@INC
);

print join "\n", @files;

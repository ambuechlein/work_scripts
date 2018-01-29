#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use File::Basename;

my $inDir = shift;
my $outDir = shift;

# $file = basename($path);
# $dir =  dirname($path);
find(\&findDB, $inDir);
sub findDB {
return unless ($_ =~ /cgb_annotation/o);
#124460252_default.old
my $dir = dirname($File::Find::name);
$dir =~ /(\d+)_default/;
return if ($dir =~ /(\d+)_default.old/);
my $pipeline_id = $1;
print "$pipeline_id\n";
my $newdir = $outDir."/$pipeline_id";
print $newdir, "\n" unless(-e $newdir);
system("mkdir -p $newdir") unless(-e $newdir);
my $stat = `cp -r $File::Find::name $newdir`; 
die "Failed to copy file:\n $stat\n" if($stat);
}

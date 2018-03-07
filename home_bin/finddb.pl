#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use File::Basename;

my $inDir = shift;

find(\&findDB, $inDir);
sub findDB {
return unless ($_ =~ /run_result/o);
my $dir = dirname($File::Find::name);
$dir =~ /(\d+)_workbench_[nuc|prot]/;
return if ($dir =~ /(\d+)_[nuc|prot]\.old/);
return if ($dir =~ /(\d+)_[nuc|prot]\.bak/);
my $pipeline_id = $1;
return unless $pipeline_id;
my $newdir = "/research/projects/isga/prod/repository/databases/$pipeline_id";
my $newdir2 = "/research/projects/isga/prod/project/output_repository/asn2all/${pipeline_id}_default/";

return if (-e $newdir2);

#print "$File::Find::name\n";
my $filename;
if($_ =~ /$pipeline_id\_run_result_nuc_db(.*)/){
$filename = $newdir.'/cgb_annotation.cds.fna' . $1;
#print $filename, "\n\n";
}else{
$_ =~ /$pipeline_id\_run_result_prot_db(.*)/;
$filename = $newdir.'/cgb_annotation.aa.fsa' . $1;
#print $filename, "\n\n";
}

print "exists: $File::Find::name\t$filename\n" if (-e $filename);
system("mkdir -p $newdir") unless(-e $newdir);
my $stat = `cp -r $File::Find::name $filename`;
die "Failed to copy file:\n $stat\n" if($stat);

}

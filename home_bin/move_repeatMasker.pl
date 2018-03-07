#! /usr/bin/perl -w
use strict;

use Getopt::Long;
my %options;
my $results = GetOptions (\%options,
                          'input_directory|i=s',
                          'prefix|p=s',
                          'output_directory|d=s',);
die "Must provide input directory\n" unless ($options{input_directory});
my $dir = $options{input_directory};
my $pre = $options{prefix};
my $out = $options{output_directory};

die "Directory does not exist\n" unless (-e $dir);
die "Directory is not a directory\n" unless (-d $dir);

opendir(my $din, "$dir") or die "Could not open directory\n";
while(my $file = readdir($din)){
  system("mv $dir/$file $out") if($file =~ /$pre/o);
#  print $cmd, "\n" if($file =~ /$pre/o);
}
closedir $din;
exit;

#! /usr/bin/perl -w
use strict;

use Getopt::Long;
use File::Find;
my %options;
my $results = GetOptions (\%options,
                          'input_directory|i=s',
                          'days|d=s');
die "Must provide input directory\n" unless ($options{input_directory});
die "Must provide days\n" unless ($options{days});
my $dir = $options{input_directory};
my $days = $options{days};

die "Directory does not exist\n" unless (-e $dir);
die "Directory is not a directory\n" unless (-d $dir);

#chdir $dir;
opendir DIR, "$dir" or die "Could not open directory\n";
foreach my $file (grep {-f && ($days < -M)} readdir DIR) {
#        unlink $file;
print $file, "\n";
}
closedir DIR;
exit;

#! /usr/bin/perl -w
use strict;

use Getopt::Long;
use File::Find;
my %options;
my $results = GetOptions (\%options,
                          'input_directory|i=s',
                         );
die "Must provide input directory\n" unless ($options{input_directory});
my $dir = $options{input_directory};

die "Directory does not exist\n" unless (-e $dir);
die "Directory is not a directory\n" unless (-d $dir);

#chdir $dir;
opendir DIR, "$dir" or die "Could not open directory\n";
foreach my $file (readdir DIR) {
  next if(not -f $file);
  next if($file =~ /\.out$/o);
#  warn $file, "\n";
  unlink $file;
}
closedir DIR;
exit;

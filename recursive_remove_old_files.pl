#! /usr/bin/perl -w
use strict;

use Getopt::Long;
use File::Find;
my %options;
my $results = GetOptions (\%options,
                          'input_directory|i=s',
                          'days|d=s',
                          'safety|s=s');
die "Must provide input directory\n" unless (defined $options{input_directory});
die "Must provide days\n" unless (defined $options{days});
my @dir;
push(@dir, $options{input_directory});
my $days = $options{days};
find(\&process_file, @dir);
exit;

sub process_file{
#    print $File::Find::name, "\n" if ($options{days} < -M $_);
  if($options{safety}){
    unlink $File::Find::name if ( ($options{days} < -M $_) and (not -d $File::Find::name) );
  } else {
    print $File::Find::name, "\n" if ($options{days} < -M $_);
  }
}

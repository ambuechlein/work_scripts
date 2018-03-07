#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "in|i=s",
                         "out|o=s",
                        );

open(my $in, $options{'in'} =~ /.gz(ip)?$/ ? "zcat $options{'in'} |" : $options{'in'} =~ /.bz(ip)?2$/ ? "bzcat $options{'in'} |" : $options{'in'}) || die("Open error: $options{'in'}");
open(my $out, ">$options{out}") or die "Can't open $options{out} for writing: $!\n";
my @querysize;
while (<$in>){
  if (/^@(\S+)/) {
         my $id = $_;
         my $seq = <$in>;
         my $id2 = <$in>;
         my $qual = <$in>;
         $seq =~ s/U/T/ig;
         print $out "$id$seq$id2$qual";
  }
}
close $in;
close $out;
exit;

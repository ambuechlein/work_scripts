#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
#use XML::Twig;
#use File::Copy;


my %options = ();
my $results = GetOptions (\%options, 'component_xml|c=s' ) || die "Error getting options?  Looks like it";

open(my $in, $options{component_xml}) or die "Can't open component.xml: $!\n";
open(my $out, ">$options{component_xml}.new") or die "Can't open component.xml.new: $!\n";
my $complete = 0;
while(<$in>){
  if($_ =~ /\<state\>/){
    $_ =~ s/\<state\>\S+\<\/state\>/\<state\>complete\<\/state\>/o;
    print $out $_;
  } elsif($_ =~ /\<total\>(\d+)\<\/total\>/) { 
    $complete = $1;
    print $out $_;
  } elsif($_ =~ /\<complete\>/) {
    $_ =~ s/\<complete>\d+\<\/complete\>/\<complete\>$complete\<\/complete\>/o;
    print $out $_;  
  } elsif($_ =~ /\<(incomplete|failed|pending|errors|running|waiting|interrupted)\>/) {
    $_ =~ s/\<(incomplete|failed|pending|errors|running|waiting|interrupted)\>\d+\<\/(incomplete|failed|pending|errors|running|waiting|interrupted)\>/\<$1\>0\<\/$2\>/;
    print $out $_;
  } else {
    print $out $_;
  } 
}
close $in;
close $out;

rename($options{component_xml}, "$options{component_xml}.orig");
rename("$options{component_xml}.new", $options{component_xml});
exit;

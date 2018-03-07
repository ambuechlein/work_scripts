#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $results = GetOptions(\%options,
                         "countstable|c=s",
                         "repeatsdir|r=s",
                        );
my $countstable = $options{countstable};;
my $repeatsdir = $options{repeatsdir};

my %map;
opendir(my $din, $repeatsdir) or die "Can't open directory $repeatsdir: $!\n";
while(my $file = readdir($din)){
  next unless($file =~ /\.nonredundant.fasta.tsv$/);
  open(my $rpt, "$repeatsdir/$file") or die "Can't open repeatMasker output $repeatsdir/$file: $!\n";
  while(<$rpt>){
    next if($_ =~ /^\s*$/o);
    chomp $_;
    $_ =~ s/^\s*//o;
    my @line = split(/\s+/, $_);
    $map{$line[9]} = $line[10];
  }
  close $rpt;
}
closedir $din;

open(my $in, $countstable) or die "Can't open countstable: $!\n";
my $h = <$in>;
my @hr = split(/\t/, $h);
shift(@hr);
print "Repeat\tRepeat Class\t", join("\t", @hr);
while(<$in>){
  my @line = split(/\t/, $_);
  my $repeat = shift(@line);
  print "$repeat\t$map{$repeat}\t", join("\t", @line);
}
close $in;
exit;

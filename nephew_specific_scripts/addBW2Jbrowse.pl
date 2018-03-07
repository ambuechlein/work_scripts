#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "sampleorder|s=s",
                        "out|o=s",
                        "jbrowse|j=s",
                       );
my $dir = $options{dir};
my $map = $options{sampleorder};
my $jb = $options{jbrowse};

open(my $out, ">$options{out}") or die "Can't open output file: $!\n";
print $out <<EOT
#!/bin/csh
EOT
;

open(my $min, $map) or die "can't open $map: $!\n";
my %labels;
my @order;
my %def;
while(<$min>){
  chomp $_;
  my @line = split(/\t/, $_);
  $line[1] =~ s/(\(|\))/_/go; $line[1] =~ s/_+/_/go;
  $labels{$line[1]} = $line[0];
  push(@order, $line[1]);
}
close $min;

print $out "\n", 'echo "Bedgraphs"', "\n";
foreach my $pre (reverse @order){
  my $nfile = "$pre.negative.bw";
  my $pfile = "$pre.positive.bw";

  warn "$nfile does not exist\n" and next unless(-e "$dir/$nfile");
  print $out "perl /nfs/labs/nephew/scripts/bw2Jbrowse.pl $jb $nfile /nfs/labs/nephew/databases/gencode_m15/chrSizes.tsv \"$pre negative\" \"$pre negative\" linear local 1 \"red\" 30 50\n"; 
  print $out "perl /nfs/labs/nephew/scripts/bw2Jbrowse.pl $jb $pfile /nfs/labs/nephew/databases/gencode_m15/chrSizes.tsv \"$pre positive\" \"$pre positive\" linear local 0 \"black\" 30 50\n";
}

close $out;
exit;

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
  my $file = "$labels{$pre}.chr.bedgraph";
  warn "$file does not exist\n" and next unless(-e "$dir/$file");
  
  print $out "perl /home/abuechle/bin/bedGraph2Jbrowse2.pl $jb $file /nfs/labs/nephew/human_databases/gencode_v26/chrSizes.tsv \"$pre positive\" \"$pre positive\" linear local 0 \"black\" 30 50\n";
}

close $out;
exit;

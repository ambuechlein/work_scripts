#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my %options;
my $results = GetOptions(\%options,
                         "desc_file|d=s",
                         "table|t=s",
                        );
my $desc_file = $options{desc_file};
my $table = $options{table};
open(my $in, "<$desc_file") or die "Can't open $desc_file.\n";

my %all_genes;
while(<$in>){
  chomp $_;
  my ($GeneID,	$desc,	$EnsFamily) = split(/\t/,$_);
  $EnsFamily = 'NA' if(not defined $EnsFamily or $EnsFamily eq '');
  $all_genes{$GeneID} = "$desc\t$EnsFamily";
}
close $in;

open(my $tb, "$table") or die "Can't open table file: $!\n";
my $header = <$tb>;
chomp $header;
my @h = split(/\t/, $header);
#     my @dwarfs = qw(Doc Grumpy Happy Sleepy Sneezy Dopey Bashful);
#    splice @dwarfs, 3, 0, 'SnowWhite';
#    print "@dwarfs";
#    # Doc Grumpy Happy SnowWhite Sleepy Sneezy Dopey Bashful
splice(@h,3,0,'Gene description','Ensembl Family Description');
my $print_h = join("\t", @h);
print "$print_h\n";

while(<$tb>){
  chomp $_;
  my @line = split(/\t/, $_);
  my $id = $line[0];
  $id =~ s/\.\d+//go;
  my $desc = "NA\tNA";
  if(defined $all_genes{$id}){
    $desc = $all_genes{$id};
  }
  splice(@line,3,0,$desc);
  my $print = join("\t", @line);
  print $print, "\n";

}

exit;


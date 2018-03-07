#!/usr/bin/env perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Aaron Buechlien
## Date: 2011

## Convert a btab blast output to tab delimeted output

use strict;
use Getopt::Long;
my ($inFile, $outFile, $helpFlag);

GetOptions(
        "h|helpFlag:s"  => \$helpFlag,
	"i|inFile:s"	=> \$inFile,		## input file
	"o|outFile:s"	=> \$outFile,		## output file
);

checkOptions();

open(my $fd, $inFile) or die "Can't open input file $inFile: $!\n";
open(my $out, ">$outFile") or die "Can't open input file $outFile: $!\n";

print $out "query_name\tdate\tquery_length\talgorithm\tdatabase_name\thit_name\tqry_start\tqry_end\thit_start\thit_end\tpercent_identity\tpercent_similarity\traw_score\tbit_score\tNULL\thit_description\tblast_frame\tqry_strand (Plus | Minus)\thit_length\te_value\tp_value\n";
#my %print_results;
<$fd>;
my $name_col=0;
my $score_col=13;
my (%max, %best, @names);
while(<$fd>){
  chomp $_;
  my @F=split /\t/, $_;
  my ($n, $s) = @F[$name_col, $score_col];
  push @names, $n if (! exists($max{$n}));
  if (! exists($max{$n}) || $s > $max{$n}) {
    $max{$n} = $s;
    $best{$n} = ();
  }
  if ($s == $max{$n}) {
    $best{$n} .= "$_\n";
  }
}
close($fd);
foreach my $k (keys %best){
  print $out $best{$k};
}
close($out);

exit;

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;

	if ($helpFlag || !defined $inFile)
	{
		die("Arguments: [-i in_file] [-o out_file]\n"
		  . "\t\n"
		  );
	}
}

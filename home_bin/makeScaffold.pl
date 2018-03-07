#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

my $usageDescription = qq(
	--input <input> \t: Input file containing the FASTA sequences.
	--output <output> \t: The output scaffold sequence file.
	--gff3 <gff3> \t\t: The GFF3 file.
	--spacer <length> \t: The number of Ns to be placed as a spacer separating the contigs.
	--verbose \t\t: Flag to output the action messages.
	--help \t\t\t: Flag to display the usage of this script.
);
my $usage = "$0 --input <input> --output <output> --gff3 <gff3> [--spacer <length>] [--verbose] [--help]\n\n$usageDescription";
my ($input, $output, $gff3, $spacerLength, $verbose, $help);
$spacerLength = 0;
my $validOptions = GetOptions ("input=s" => \$input, "output=s" => \$output, "gff3=s" => \$gff3, "spacer=i" => \$spacerLength, "verbose"  => \$verbose, "help" => \$help);
if ($help) {
	print "$usage\n";
	exit 0;
}

unless ($validOptions &&  $input && $output && $gff3 && $spacerLength) {
	print "SYNTAX ERROR\n$usage\n";
	exit 1;
}

if ( $verbose ) {
	print "Checking Permissions...\n";
}

unless (open INPUT, '<', $input) {
	print "ERROR: Unable to open the file $input.\n";
	exit 1;
}

unless (open OUTPUT, '>', $output) {
	print "ERROR: Unable to create the file $output.\n";
	exit 1;
}

print OUTPUT ">SeqID\n";

unless (open GFF3, '>', $gff3) {
	print "ERROR: Unable to create the file $gff3.\n";
	exit 1;
}

print GFF3 "##gff-version 3\n";

if ($verbose) {
	print "done.\n";
}

if ( $verbose ) {
	print "Loading sequences from $input... \n";
}

my %sequences = ();
my %lengths = ();
my %annotations = ();
my ($name, $sequence, $annotation) = ('', '', '');

while (<INPUT>) {
	chomp;
	$_ = trim($_);
	if ($_ ne '')  {
		if ((substr $_, 0, 1) eq ">") {
			if ( $name && $sequence ) {
				$sequences{$name} = $sequence;
				$lengths{$name} = length $sequence;
				$annotations{$name} = $annotation;
				$sequence = '';
				$name = '';
				$annotation = '';
			}
			$_ =~ s/>//;
			my @headerElements = split /\s+/;
			$name = $headerElements[0];

			$annotation="ID=$name;Name=$name";
			for (my $i = 1; $i < scalar @headerElements; $i++) {
				$annotation .= ";$headerElements[$i]";
			}
		} else {
			$sequence .= $_;
		}
	}
}

if ( $name && $sequence ) {
	$sequences{$name} = $sequence;
	$lengths{$name} = length $sequence;
	$annotations{$name} = $annotation;
	$sequence = '';
	$name = '';
	$annotation = '';
}
close INPUT;

if ( $verbose ) {
	print "Finished loading sequences from $input.\n";
}

my $spacer = 'N' x $spacerLength;

if ($verbose) {
	print "Sorting sequences on lengths... \n";
}

sub hashValueDescendingNum {
	$lengths{$b} <=> $lengths{$a};
}

my $start = 1;

foreach my $key (sort hashValueDescendingNum (keys(%lengths))) {
	print OUTPUT $sequences{$key};
	print GFF3 "SeqID\tmakeScaffold\tcontig\t$start\t";
	print GFF3 ($start + $lengths{$key} - 1);
	$start += $lengths{$key};
	print GFF3 "\t.\t.\t.\t$annotations{$key}\n";
	print OUTPUT $spacer;
	print GFF3 "SeqID\tmakeScaffold\tspacer\t$start\t";
	print GFF3 ($start + $spacerLength - 1);
	$start += $spacerLength;
	print GFF3 "\t.\t.\t.\t.\n";
}
print OUTPUT "\n";

if ($verbose) {
	print "Done. \n";
}
print GFF3 "SeqID\tmakeScaffold\tscaffold\t1\t";
print GFF3 $start - 1;
print GFF3 "\t.\t.\t.\tID=SeqID;Name=SeqID\n";
close GFF3;

if ($verbose) {
	print "Generating output... \n";
}
close OUTPUT;

if ($verbose) {
	print "done.\n";
}

exit;


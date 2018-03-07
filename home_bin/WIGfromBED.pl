#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use List::Util 'max';
use File::Basename;

sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

my $usageDescription = qq(
	--name \t\t: The name part of the WIG header
	--description\t: The description part of the WIG header
	--bedFile\t: The input bed file
	--wigFile\t: The output wiggle file
	--verbose \t: Flag to output the action messages.
	--help \t\t: Flag to display the usage of this script.
);

my $usage = "$0 --bedFile <bedFile> --wigFile <wigFile> [--name <WIG Name>] [--description <WIG Description>] [--verbose] [--help]\n\n$usageDescription";

my ($bedFile, $wigFile, $name, $description, $verbose, $help);

my $validOptions = GetOptions ("bedFile=s" => \$bedFile, "wigFile=s" => \$wigFile, "name=s" => \$name, "description=s" => \$description, "verbose"  => \$verbose, "help" => \$help);

if ($help) {
	print "$usage\n";
	exit 0;
}

unless ($validOptions &&  $bedFile && $wigFile) {
	print "SYNTAX ERROR\n$usage\n";
	exit 1;
}

$name = basename($bedFile) unless ($name);
$description = $name unless ($description);


if ( $verbose ) {
	print "Reading $bedFile\n";
}
	
unless (open(BED, '<', $bedFile)) {
       	print "ERROR: Unable to open the input file $bedFile.\n";
        exit 1;
}

unless (open(WIG, '>', $wigFile)) {
       	print "ERROR: Unable to open the output file $wigFile.\n";
        exit 1;
}

print WIG 'track type=wiggle_0 name="' . $name . '" description="' . $description . '"' . "\n";
my %data = ();

while (<BED>) {
	chomp;
	if ((trim($_) ne '') && ((substr $_, 0, 1) ne "#")) {
		my @bedData = split /\t+/;
		my ($chr, $start, $end) = ($bedData[0], $bedData[1] + 1, $bedData[2]);
		for (my $i = $start; $i <= $end; $i++) {
			$data{$chr}{$i} = 0 unless ($data{$chr}{$i});
			$data{$chr}{$i}++;
		}

	}
}
close BED;

foreach my $chr(sort keys %data) {
	if ( $verbose ) {
		print "Generating wiggle data for $chr...\n";
	}
	print WIG "fixedStep chrom=$chr start=1 step=1 span=1\n";
	my $chrLen = max(keys %{$data{$chr}});
	for (my $i = 1; $i <= $chrLen; $i++) {
		if ($data{$chr}{$i}) {
			print WIG "$data{$chr}{$i}\n";
		} else {
			print WIG "0\n";
		}
	}
}

if ( $verbose ) {
	print "Creating the file $wigFile\n";
}
close WIG;

if ( $verbose ) {
	print "Done.\n";
}


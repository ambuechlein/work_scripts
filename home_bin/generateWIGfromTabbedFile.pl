#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Getopt::Long;

sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

my $usageDescription = qq(
	--inputFile <inputFile> \t: The input tab-delimited file.
	--outputFile <outputFile> \t: The output WIG file, defaults to the extension .wig appended to the inputFile name.
	--chr <int> \t\t\t: Column representing the chromosome.
	--pos <int> \t\t\t: Column representing the position, not to be used along with startPos and endPos options.
	--begin <int> \t\t\t: Column representing the start position.
	--end <int> \t\t\t: Column representing the end position.
	--BedStart \t\t\t: Actual start one position after begin, as in BED format.
	--score <int> \t\t\t: Column representing the score.
	--name <name> \t\t\t: Name to be used in the header, defaults to the name of the inputFile.
	--description <description> \t: Description to be used in the header, defaults to the name.
	--verbose \t\t\t: Flag to output the action messages.
	--help \t\t\t\t: Flag to display the usage of this script.

	Note: Columns are zero based, start with 0 and not 1.
);
my $usage = "$0 --inputFile <inputFile> [--outputFile <outputFile>] --chr <int> [--pos <int>] [--begin <int> --end <int>] [--BedStart] --score <score> [--name <name>] [--description <description>]  [--verbose] [--help]\n\n$usageDescription";

my ($inputFile, $outputFile, $colChr, $colPos, $colBegin, $colEnd, $bedStart, $colScore, $name, $description, $verbose, $help);

my $validOptions = GetOptions ("inputFile=s" => \$inputFile, "outputFile=s" => \$outputFile, "chr=i" => \$colChr, "pos=i" => \$colPos, "begin=i" => \$colBegin, "end=i" => \$colEnd, "BedStart" => \$bedStart, "score=i" => \$colScore, "name=s" => \$name, "description=s" => \$description, "verbose"  => \$verbose, "help" => \$help);

if ($help) {
	print "$usage\n";
	exit 0;
}

unless ($validOptions && defined $inputFile && defined $colChr && $colChr >= 0 && ((defined $colPos && $colPos >= 0) || (defined $colBegin && defined $colEnd && $colBegin >= 0 && $colEnd >= 0)) && defined $colScore && $colScore >= 0) {
#unless ($validOptions && $inputFile && $colChr && ($colPos || ($colBegin && $colEnd)) && $colScore) {
	print "SYNTAX ERROR\n$usage\n";
	exit 1;
}

$outputFile = "$inputFile.wig" unless ($outputFile);
$name = basename($inputFile) unless ($name);
$description = $name unless ($description);

print "Checking for file permissions... " if ($verbose);

unless (open(INPUT, '<', $inputFile)) {
	print "\nERROR: Unable to open the file $inputFile.\n";
	exit 1;
}

my $openLine = <INPUT>;
chomp $openLine;
       	
unless ($openLine =~ /^track type=wiggle_0/) {
	$openLine = 'track type=wiggle_0 name="' . $name . '" description="' . $description . '"';
	close INPUT;
	unless (open(INPUT, '<', $inputFile)) {
		print "\nERROR: Unable to open the file $inputFile.\n";
		exit 1;
	}
}

unless (open(OUTPUT, '>', $outputFile)) {
	print "\nERROR: Unable to create the file $outputFile.\n";
	exit 1;
}

print "Done!\n" if ($verbose);
	
print OUTPUT "$openLine\n";
my $oldChr = '';
my $oldEnd = 0;

print "Reading data from $inputFile...\n" if ($verbose);

while (<INPUT>) {
	chomp;
	if (($_ ne '') && ((substr $_, 0, 1) ne "#")) {
		my @inputDataCols = split /\t+/;
		my ($chr, $start, $end, $score);
		$chr = $inputDataCols[$colChr];
		if ($colPos) {
			$start = $inputDataCols[$colPos];
			$end = $inputDataCols[$colPos];
		} else {
			$start = $inputDataCols[$colBegin];
			$end = $inputDataCols[$colEnd];
			$start++ if ($bedStart);
		}
		$score = $inputDataCols[$colScore];
		if ($oldChr ne $chr) {
			print "Generating wiggle output for $chr\n" if ($verbose);
			print OUTPUT "fixedStep chrom=$chr start=1 step=1 span=1\n";
			$oldChr = $chr;
			$oldEnd = 0;
		}
		unless ($start == $oldEnd + 1) {
			for (my $j = $oldEnd + 1; $j < $start; $j++) {
				print OUTPUT "0\n";
			}
		}
		for (my $i = $start; $i <= $end; $i++) {
			print OUTPUT "$score\n";
		}
		$oldEnd = $end;
	}
}

close OUTPUT;
close INPUT;

print "$outputFile created.\n" if ($verbose);



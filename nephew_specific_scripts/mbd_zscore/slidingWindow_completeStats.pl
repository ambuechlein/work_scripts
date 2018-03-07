#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil);
use Getopt::Long;

sub mean {
	my($data) = [@_]; 
	die("Empty array\n") if (scalar @$data == 0); 
	my $total = 0; 
	foreach (@$data) { 
		$total += $_; 
	} 
	my $average = $total / @$data; 
	return $average;
}
sub standard_deviation{
	my($data) = [@_]; 
	return 0 if(scalar @$data == 1);
	my $average = &mean(@$data); 
	my $sqtotal = 0; 
	foreach(@$data) { 
		$sqtotal += ($average-$_) ** 2; 
	} 
	my $std = ($sqtotal / (@$data-1)) ** 0.5; 
	return $std; 
}

my $usageDescription = qq(
	--sizes <file>\t: Chromosome sizes file
	--input <file>\t: Input data file, default STDIN
	--chr_col <int>\t: Chromosome column in the data file
	--pos_col <int>\t: Position column in the data file
	--start_col <int>\t: Start position column in the data file
	--end_col <int>\t: End position column in the data file
	--BedStart \t: Start positions are zero based, as in BED format files
	--score_col <int>\t: Score column in the data file
	--output <file>\t: The output file with window scores, default STDOUT
	--winsize <int>\t: Size of the window, default 100
	--stepsize <int>\t: Step size for the slide, default 20
	--verbose \t: Flag to output action messages
	--help \t\t: Flag to display the usage of this script

	The column numbers start from 1 and not 0
);

my $usage = "$0 --sizes <file> [--input <file>] --chr_col <int> (--pos_col <int> | --start_col <int> --end_col <int> [--BedStart]) --score_col <int> [--output <file>] [--verbose] [--help]\n\n$usageDescription";

my ($sizes, $input, $chr_col, $pos_col, $start_col, $end_col, $bedStart, $score_col, $output, $winsize, $stepsize, $verbose, $help);

my $validOptions = GetOptions ("sizes=s" => \$sizes, "input=s" => \$input, "chr_col=i" => \$chr_col, "pos_col=i" => \$pos_col, "start_col=i" => \$start_col, "end_col=i" => \$end_col, "BedStart" => \$bedStart, "score_col=i" => \$score_col, "output=s" => \$output, "winsize=i" => \$winsize, "stepsize=i" => \$stepsize, "verbose"  => \$verbose, "help" => \$help);


if ($help) {
	print STDERR "$usage\n";
	exit 0;
}

unless ($validOptions &&  defined $sizes && $chr_col && ($pos_col || ($start_col && $end_col)) && $score_col) {
	print STDERR "SYNTAX ERROR\n$usage\n";
	exit 1;
}
	
print STDERR "Checking for file permissions...\n" if (defined $verbose);

unless (open(SIZES, '<', $sizes)) {
	print STDERR "ERROR: Unable to open the file $sizes.\n";
	exit 2;
}

if (defined $input) {
	if ($input =~ /\.gz$/) {
		unless (open INPUT, "gunzip -c $input |") {
			print STDERR "ERROR: Unable to open the file $input.\n";
			exit 3;
		}
	} else {
		unless (open(INPUT, '<', $input)) {
			print STDERR "ERROR: Unable to open the file $input.\n";
			exit 4;
		}
	}
} else {
	*INPUT = *STDIN;
}

if ($output) {
	unless (open(OUTPUT, '>', $output)) {
		print STDERR "ERROR: Unable to create the file $output.\n";
		exit 5;
	}
} else {
	*OUTPUT = *STDOUT;
}
print OUTPUT "#Chr\tBin_ID\tBin_Start\tBin_End\tAverage_Bin_Score\tStandard_Deviation\tIndex\n";

$chr_col--;
if (defined $pos_col) {
	$pos_col--;
} else {
	$start_col--;
	$end_col--;
}
$score_col--;

$winsize = 100 unless($winsize);
$stepsize = 20 unless($stepsize);

print STDERR "Done.\n" if (defined $verbose);

my %chrSizes = ();

print STDERR "Loading chromosome sizes...\n" if (defined $verbose);
my $size_counter = 0;
while(<SIZES>) {
	chomp;
	next if (/^#/);
	my @cols = split /\t/;
	$chrSizes{$cols[0]} = $cols[1];
	$size_counter++;
	print STDERR "Loaded $size_counter chromosomes.\r" if ($verbose);
}
close SIZES;
print STDERR "Total $size_counter chromosomes.            \n\n" if (defined $verbose);

my %data = ();

print STDERR "Loading data... \n" if (defined $verbose);
my $data_counter = 0;
while(<INPUT>){
	chomp;
	next if (/^#/);
	my @cols = split /\t/;
	$data{$cols[$chr_col]}{$cols[$pos_col]} = $cols[$score_col];
	$data_counter++;
	if ($data_counter%1000 == 0) {
		print STDERR "Loaded $data_counter rows.\r" if ($verbose);
	}
}
close INPUT;
print STDERR "Total $data_counter rows.                \n\n" if (defined $verbose);

my $mid = ceil($winsize / 2);

print STDERR "Computing bin averages... \n\n" if (defined $verbose);
foreach my $chr (sort keys %chrSizes) {
	my @activememships = ();
	my $nextbin = $mid;
	my %binscores = ();
	print STDERR "Starting with $chr of size $chrSizes{$chr} bases.\n\n" if (defined $verbose);
	for my $n (1..$chrSizes{$chr}){
		if ($n + $mid > $nextbin) {
			if ($nextbin + $mid <= $chrSizes{$chr}) {
				push @activememships, $nextbin;
				$nextbin += $stepsize;
			}
		}
		if (scalar @activememships > 0) {
			shift @activememships if (($n - $activememships[0]) > $mid);
			if (scalar @activememships > 0) {
				foreach my $memshipid (@activememships) {
					$binscores{$memshipid} = () unless (exists $binscores{$memshipid});
					my $score = 0;
					$score = $data{$chr}{$n} if (exists $data{$chr}{$n});
					push @{$binscores{$memshipid}}, $score;
				}
			}
		}
		if ($n%1000 == 0) {
			print STDERR "Scanned $n bases.\r" if ($verbose);
		}
	}
	print STDERR "Scanned all the bases in chromosome $chr.\n" if ($verbose);
	print STDERR "Generating bin statistics for chromosome $chr.\n" if ($verbose);
	my $index = 0;
	foreach my $bin_id (sort {$a<=>$b} keys %binscores) {
		$index++;
		my $average = 0;
		my $std = 0;
		$average = &mean(@{$binscores{$bin_id}});
		$std = &standard_deviation(@{$binscores{$bin_id}});
		print OUTPUT "$chr\t$bin_id\t" . ($bin_id - $mid + 1) . "\t" . ($bin_id + $mid) . "\t$average\t$std\t$chr" . "_$index\n";
	}
	print STDERR "\n\n" if ($verbose);
}
close OUTPUT;



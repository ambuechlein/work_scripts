#!/usr/bin/env perl

my $fastq = shift;
my $fivePrimeAdapter = shift;
my $max = shift;
my $min = shift;
my $min_length = shift;
my %stats;
$max = 10 unless defined($max);
$min = 1 unless defined($min);
$min_length = 17 unless defined($min_length);
my $cut = 0;
if (length($fivePrimeAdapter) < 6) {
	print STDERR "ERROR: five prime adapter \"$fivePrimeAdapter\" is too short. Must be $max bases or longer.\n";
	exit 100;
}

my @adapt = split //,$fivePrimeAdapter;
open(my $ffh, $fastq) or die "Can't open $fastq: $!\n";
while (<$ffh>) {
	chomp;
	my $id = $_;
	my $seq = <$ffh>;
	my $sep = <$ffh>;
	my $qual = <$ffh>;
	chomp $seq;
	chomp $qual;
#	print "$id\n$seq\n$sep$qual\n";
	my $len = length($seq);
	my @seq = split //,$seq;
	my $match;
	for (my $i = 0; $i < ($max-$min+1); $i++) {
		my $match = 0;
		my $mismatch = 0;
		$good = 0;
		for (my $j = $i; $j < $max; $j++) {
			my $pos = $len - $max + $j;
			if ($seq[$pos] eq $adapt[$j-$i]) {
				$match++;
			} else {
				$mismatch++;
				last if $mismatch > 1;
			}
#			print "$i\t$j\t$match\t$pos\t$seq[$pos]\t$adapt[$j-$i]\n";
		}
#		my $v = $max - $i - 1;
#		print "$i\t$max\t$match\t$v\t$len\n";
		if (($max-$i) >= 6) {
			$good = 1 if $match >= ($max - $i - 1);
		} else {
			$good = 1 if $match == ($max - $i);
		}
		if ($good) {
			$seq = join "",@seq[0..($len-($max-$i+1))];
			my @qual = split //,$qual;
			$qual = join "",@qual[0..($len-($max-$i+1))];
			if (length($seq) >= $min_length) {
				print "$id\n$seq\n$sep$qual\n";
				$stats{$max-$i}++;
			} else {
				$cut++;
			}
			last;
		}
	}
	if ($good) {	
	} else {
		print "$id\n$seq\n$sep$qual\n";
		$stats{0}++;
	}
#	print "$match\t$good\n";
#	print "$id\n$seq\n$sep$qual\n";
}

for (my $i = 0; $i <= $max; $i++) {
	my $v = $stats{$i}+0;
	print STDERR "$i\t$v\n";
}
print STDERR "cut\t$cut\n";

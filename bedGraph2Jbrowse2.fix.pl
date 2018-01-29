#!/usr/bin/env perl
use strict;
use warnings;
use FileUtil;

#perl ~/bin/bedGraph2Jbrowse2.pl /nfs/bio/web/tools/jbrowse/databases/mbd_07-15-14/ /nfs/labs/nephew/a2780_methylation/mbd_07-15-14/rita/110829_s_1_sorted.txt.hg19.sorted.bedgraph /nfs/bio/db/Homo_sapien/hg19_all/chrSizes.tsv "110829_s_1_sorted" "110829_s_1_sorted" linear local 0 "black" 30 50
my $jbrowseDir = shift;
my $bedGraphFile = shift;
my $sizeFile = shift;
my $label = shift;
my $key = shift;
my $scale = shift;
my $scope = shift;
my $minus = shift;
my $posColor = shift;
my $height = shift;
my $linearMax = shift;
my $negColor = $posColor;

my $bwFile = $bedGraphFile;
$bwFile =~ s/.*\///;
$bwFile =~ s/bedgraph$/bw/i;

$posColor = "black" unless defined($posColor);
$negColor = "red" unless defined($negColor);

if ($height) {
} else {
	if ($scale eq "log") {
		$height = 50;
	} else {
		$height = 100;
	}
}

my $cmd;

my @x;
for (my $i = 6; $i < 12; $i++) {
	$x[$i] = " " x $i;
}

if (-e $bedGraphFile && -s $bedGraphFile) {
} else {
	print STDERR "ERROR: required file cannot be found or is empty: $bedGraphFile\n";
	exit 1;
}
if (-e $sizeFile && -s $sizeFile) {
} else {
	print STDERR "ERROR: required file cannot be found or is empty: $sizeFile\n";
	exit 3;
}

#$cmd = "mkdir -p $jbrowseDir/bigWig";
#runSystem($cmd);
my $logMax = 10;
if ($minus) {
	my $tempFile = "/tmp/tmp.$$.negBed";
	my $ofh = new FileHandle(">$tempFile");
	my $tfh = FileUtil::openFileHandle($bedGraphFile);
	while (<$tfh>) {
		if (/^track/i) {
			print $ofh $_;
		} else {
			chomp;
			my @d = split /\t/;
			$logMax = $d[3] if $d[3] > $logMax;
			$d[3] = -$d[3];
			my $l = join "\t",@d;
			print $ofh $l,"\n";
		}
	}
	$ofh->close();
	$tfh->close();
	$cmd = "/home/abuechle/bin/bedGraphToBigWig $tempFile $sizeFile $jbrowseDir/bigWig/$bwFile";
	runSystem($cmd);
	sleep 1;
	$cmd = "rm $tempFile";
	runSystem($cmd);
} else {
	my $tempFile = "/tmp/tmp.$$.posBed";
	my $ofh = new FileHandle(">$tempFile");
	my $tfh = FileUtil::openFileHandle($bedGraphFile);
	while (<$tfh>) {
		if (/^track/i) {
			print $ofh $_;
		} else {
			chomp;
			my @d = split /\t/;
			$logMax = $d[3] if $d[3] > $logMax;
			my $l = join "\t",@d;
			print $ofh $l,"\n";
		}
	}
	$ofh->close();
	$tfh->close();
	$cmd = "/home/abuechle/bin/bedGraphToBigWig $tempFile $sizeFile $jbrowseDir/bigWig/$bwFile";
	runSystem($cmd);
	sleep 1;
	$cmd = "rm $tempFile";
	runSystem($cmd);
}

sub runSystem {
  my ($cmd,$message) = @_;

  my $error;
  print STDERR $cmd,"\n";
  if ($error = system($cmd)) {
	  print STDERR "ERROR: $cmd\n";
    exit 1;
  }
}

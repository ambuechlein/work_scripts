use FileUtil;

my $jbrowseDir = shift;
my $bwFile = shift;
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

$posColor = "black" unless defined($posColor);
$negColor = "red" unless defined($negColor);

unless($height) {
	if ($scale eq "log") {
		$height = 50;
	} else {
		$height = 100;
	}
}

my $originalTrackFile = "$jbrowseDir/trackList.json.orig";
my $oldTrackFile = "$jbrowseDir/trackList.json.old";
my $currentTrackFile = "$jbrowseDir/trackList.json";
my $cmd;

my @x;
for (my $i = 6; $i < 12; $i++) {
	$x[$i] = " " x $i;
}

unless(-e $bwFile && -s $bwFile) {
	print STDERR "ERROR: required file cannot be found or is empty: $bedGraphFile\n";
	exit 1;
}
unless(-e $currentTrackFile && -s $currentTrackFile) {
	print STDERR "ERROR: required file cannot be found or is empty: $currentTrackFile\n";
	exit 2;
}
unless(-e $sizeFile && -s $sizeFile) {
	print STDERR "ERROR: required file cannot be found or is empty: $sizeFile\n";
	exit 3;
}

$cmd = "mkdir -p $jbrowseDir/bigWig";
runSystem($cmd);
$cmd = "cp $bwFile $jbrowseDir/bigWig/";
runSystem($cmd);

unless (-e $originalTrackFile && -S $originalTrackFile) {
	$cmd = "mv $currentTrackFile $originalTrackFile";
	runSystem($cmd);
	$cmd = "cp $originalTrackFile $oldTrackFile";
	runSystem($cmd);
} else {
	$cmd = "mv $currentTrackFile $oldTrackFile";
	runSystem($cmd);
}

my $qq = chr(34);
my $ofh = new FileHandle(">$currentTrackFile");
my $fh = FileUtil::openFileHandle($oldTrackFile);
$logMax = int($logMax/10+1) * 10;
while (<$fh>) {
	if (/^\s*\"tracks\"\s*:\s*\[/o) {
		print $ofh $_;
		print $ofh "$x[6]\{\n";
		print $ofh "$x[8]$qq","label","$qq : $qq$label$qq,\n";
		print $ofh "$x[8]$qq","key","$qq : $qq$key$qq,\n";
		print $ofh "$x[8]$qq","storeClass","$qq : $qq","JBrowse/Store/SeqFeature/BigWig$qq,\n";
		print $ofh "$x[8]$qq","urlTemplate","$qq : $qq","bigWig/$bwFile$qq,\n";
		print $ofh "$x[8]$qq","type","$qq : $qq","JBrowse/View/Track/Wiggle/XYPlot$qq,\n";
		print $ofh "$x[8]$qq","variance_band","$qq : false,\n";
		print $ofh "$x[8]$qq","scale","$qq : $qq","$scale$qq,\n";
		print $ofh "$x[8]$qq","autoscale","$qq : $qq","$scope$qq,\n";
		if ($scope eq "local") {
		} else {
			if ($scale eq "log") {
				if ($minus) {
					print $ofh "$x[8]$qq","min_score","$qq : $qq","-$logMax$qq,\n";
					print $ofh "$x[8]$qq","max_score","$qq : $qq","0$qq,\n";
				} else {
					print $ofh "$x[8]$qq","min_score","$qq : $qq","0$qq,\n";
					print $ofh "$x[8]$qq","max_score","$qq : $qq","$logMax$qq,\n";
				}
			}
			if ($linearMax) {
				if ($scale eq "linear") {
					if ($minus) {
						print $ofh "$x[8]$qq","min_score","$qq : $qq","-$linearMax$qq,\n";
						print $ofh "$x[8]$qq","max_score","$qq : $qq","0$qq,\n";
					} else {
						print $ofh "$x[8]$qq","min_score","$qq : $qq","0$qq,\n";
						print $ofh "$x[8]$qq","max_score","$qq : $qq","$linearMax$qq,\n";
					}
				}
			}
		}
		print $ofh "$x[8]$qq","style","$qq : {\n";
		print $ofh "$x[10]$qq","height","$qq : $height,\n";
		print $ofh "$x[10]$qq","pos_color","$qq : $qq$posColor$qq,\n";
		print $ofh "$x[10]$qq","neg_color","$qq : $qq$negColor$qq,\n";
		print $ofh "$x[8]}\n";
		print $ofh "$x[6]},\n";
	} else {
		print $ofh $_;
	}
}
$ofh->close();
$fh->close();


sub runSystem {
  my ($cmd,$message) = @_;

  my $error;
  print STDERR $cmd,"\n";
  if ($error = system($cmd)) {
	  print STDERR "ERROR: $cmd\n";
    exit 1;
  }
}

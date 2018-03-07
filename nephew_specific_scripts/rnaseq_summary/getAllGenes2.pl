use FileUtil;

my $f = shift;
my $fo = shift;
my $fh = FileUtil::openFileHandle($f);
my $fho = new FileHandle(">$fo");
while (<$fh>) {
#        next unless /TYPE:1/;
	next unless /XT:i:1/;
	next if /XL\:Z\:tRNA/;
	if (/unique_match.(\S+)/) {
#		$genes{$1}++;
#		print "$1\n";
		my $match = $1;
		my @set = split /==/,$match;
		my %good;
		foreach my $i (@set) {
			my @d = split /~~/,$i;
#print "AAA $i\t$d[1]\n";
			next unless $d[1];
			if ($i =~ /(\S+);/) {
#print "RRR $i $1\n";
				$good{$1}++;
			}
		}
#		my $x = join " ",@good;
#		print "bbb ",$x,"\n";
		foreach my $i (keys %good) {
			$genes{$i}++;
#			print "xxx $i\t$d[0]\n";
		}
	}
	if (/ambiguous_match.(\S+)/) {
		my $match = $1;
		my @set = split /==/,$match;
		my %good;
		foreach my $i (@set) {
#print $i,"\n";
			my @d = split /~~/,$i;
			next unless $d[1];
			if ($i =~ /(\S+);/) {
#print $i,"\n";
				$good{$1}++;
			}
		}
		foreach my $i (keys %good) {
			$genes{$i}++;
#			print "$i\t$d[0]\n";
		}
	}
}
$fh->close();

foreach my $i (keys %genes) {
	print $fho "$i\t$genes{$i}\n";
#	print "$i\t$genes{$i}\n";
}
$fho->close();

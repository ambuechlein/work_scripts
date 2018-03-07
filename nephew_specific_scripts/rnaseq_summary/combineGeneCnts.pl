use FileUtil;

my $cntList = shift;
my @cnts = split /,/,$cntList;

my @headers;
push @headers,"gene";
foreach my $f (@cnts) {
	my $fh = FileUtil::openFileHandle($f);
#	$f =~ /(\S+).hist$/;
#	$f =~ /(\S+).gene\.cnts$/;
#	$f =~ /(\S+).gene2\.cnts$/;
#	$f =~ /(\S+).stranded.cnts$/;
#	$f =~ /(\w+).\S+\.cnts$/;
#	$f =~ /.*\.(\S+)\.cnts$/;
#	$f =~ /(\S+)\.cnts$/;
	$f =~ /(\S+)$/;
#	$f =~ /(\w+?).stranded.cnts/;
	my $id = $1;
	push @headers,$id;
	while (<$fh>) {
		chomp;
		my @d = split /\s+/;
		$map{$d[0]}{$id} = $d[1];
#		$map{$d[0]}{$id} = $d[1];
	}
	$fh->close();
}

my $l = join "\t",@headers;
print $l,"\n";
foreach my $i (keys %map) {
	my @l;
	push @l,$i;
	for (my $j = 1; $j < @headers; $j++) {
		push @l,int(0.5+$map{$i}{$headers[$j]});
	}
	my $l = join "\t",@l;
	print $l,"\n";
}

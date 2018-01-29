use FileUtil;

my $data = shift;
my $dfh = FileUtil::openFileHandle($data);
while (<$dfh>) {
	chomp;
	my ($id,$bam,$cnts,$rm,$genes,$genome,$trans,$unmapped,$out,$tophatCmd,$screen) = split /\t/;
        warn "ID $id\nBam $bam\nCNTS $cnts\nRM $rm\nGenes $genes\nGenome $genome\nTrans $trans\nUnmapped $unmapped\nout $out\ntophatCMD $tophatCmd\nscreen $screen\n";

	my %repeat;
	my %rclass;
	if (-e $rm) {
		print STDERR "$id: Starting Repeats $rm\n";
		my $rfh = FileUtil::openFileHandle($rm);
		while (<$rfh>) {
			s/^\s+//;
			chomp;
			my @d = split /\s+/;
			$repeat{$d[4]} = $d[9];
			$rclass{$d[4]} = $d[10];
		}
		$rfh->close();
	} else {
		print STDERR "ERROR: could not find repeat masker file for $id: $rm\n";
		next;		
	}
	if (0) {
		if (-e $screen) {
			print STDERR "$id: Starting Screen $screen\n";
			my $sfh = FileUtil::openFileHandle($screen);
			while (<$sfh>) {
				chomp;
				my @d = split /\t/;
				$repeat{$d[0]} = "Ribo" unless exists($repeat{$d[0]});
				$rclass{$d[0]} = "Ribo" unless exists($rclass{$d[0]});
			}
			$sfh->close();
		} else {
			print STDERR "ERROR: could not find screen file for $id: $screen\n";
			next;		
		}
	}
	my %cnts;
	my $totalCnt;
	if (-e $cnts) {
		print STDERR "$id: Starting Counts $cnts\n";
		my $cfh = FileUtil::openFileHandle($cnts);
		while (<$cfh>) {
			chomp;
			my @d = split /\t/;
			$cnts{$d[1]}{$d[0]}=0;
			$totalCnt++;
		}
		$cfh->close();
	} else {
		print STDERR "ERROR: could not find counts file for $id: $cnts\n";
		next;		
	}
	my %types;
	my %genes;
	if (-e $genes) {
		print STDERR "$id: Starting Genes $genes\n";
		my $gfh = FileUtil::openFileHandle($genes);
		while (<$gfh>) {
#		print $_;
			chomp;
			my @d = split /\t/;
			my $type;
			if ($d[1] =~ /^\d+$/ or $d[1] =~ /^JH/ or $d[1] =~ /^GL/ or $d[1] eq 'X' or $d[1] eq 'Y' or $d[1] eq 'MT') {
				$type = 1;	### human
			} else {
				$type = 2;	### bacterial
			}
			$types{$d[0]} = $types{$d[0]} | $type;
			if ($d[3] eq "Gene") {
				if ($d[4] eq "multi_match") {
					$genes{$d[0]} = \@d;
				}
			} else {
				$genes{$d[0]} = \@d;
			}
		}
		$gfh->close();
	} else {
		print STDERR "ERROR: could not find genes file for $id: $genes\n";
		next;		
	}
	if (0) {
		my %gbest;
		my %gbact;
		my %lens;
		if (-e $genome) {
			my $gfh = FileUtil::openFileHandle($genome);
			while (<$gfh>) {
				chomp;
				my @d = split /\t/;
				/NM:i:(\d+)/;
				my $edit = $1;
				$lens{$d[0]} = length($d[9]);
				if (exists($gbest{$d[0]})) {
					$gbact{$d[0]}++;
				} else {
					$gbest{$d[0]} = $edit;
					$gbact{$d[0]}++;
				}
			}
			$gfh->close();
		} else {
			print STDERR "ERROR: could not find genome file for $id: $genome\n";
			next;		
		}
	}
	my %seen;
	my $fo = new FileHandle("| gzip -c > $out");
	if (-e $bam) {
		print STDERR "$id: Starting Hits $bam\n";
		my $bfh = FileUtil::openFileHandle("samtools view $bam|");
		while (<$bfh>) {
			chomp;
			my @d = split /\t/;
			my $id = shift @d;
			next if exists($seen{$id});
			if ($d[1] eq "chrR") {
				if (exists($repeat{$id})) {
				} else {
					$repeat{$id} = "Ribo";
					$rclass{$id} = "Ribo";
				}
			}
			/NM:i:(\d+)/;
			my $edit = $1;
			push @d,"SEQ:$id";
			my $num = keys %{$cnts{$id}};
			if (exists($rclass{$id})) {
				push @d,"REP:$repeat{$id}","RCL:$rclass{$id}";
			}
			my $einfo;
			if (exists($genes{$id})) {
				$einfo = join "!",@{$genes{$id}};
				push @d,"GENE:$einfo";
			}
			if ($einfo =~ /multi_match/) {
				$seen{$id} = 1; 
			}
			my $v = $types{$id} + 0;
			push @d,"TYPE:$v";
			foreach my $i (keys %{$cnts{$id}}) {
				my $l = join "\t",$i,@d;
#				print $l,"\n";
				print $fo $l,"\n";
			}
		}
		$bfh->close();
	} else {
		print STDERR "ERROR: could not find bam file for $id: $bam\n";
		next;		
	}
	if (-e $unmapped) {
		print STDERR "$id: Starting Unmapped $unmapped\n";
		my $ufh = FileUtil::openFileHandle("samtools view $unmapped|");
		while (<$ufh>) {
			chomp;
			my @d = split /\t/;
			my $id = shift @d;
			/NM:i:(\d+)/;
			my $edit = $1;
			push @d,"SEQ:$id";
			my $num = keys %{$cnts{$id}};
			if (exists($rclass{$id})) {
				push @d,"REP:$repeat{$id}","RCL:$rclass{$id}";
			}
			if (exists($genes{$id})) {
				my $einfo = join "!",@{$genes{$id}};
				push @d,"GENE:$einfo";
			}
			my $v = $types{$id} + 0;
			push @d,"TYPE:$v";
			foreach my $i (keys %{$cnts{$id}}) {
				my $l = join "\t",$i,@d;
				print $fo $l,"\n";
			}
		}
		$ufh->close();
	} else {
		print STDERR "ERROR: could not find bam file for $id: $bam\n";
		next;		
	}
	$fo->close();
	
}
$dfh->close();

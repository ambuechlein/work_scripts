use FileUtil;

my $list = shift;
my $order = shift;
my $gtf = shift;
my $summaryOut = shift;
my $baseHistOut = shift;

### Globals
my %totalReads;
my %lengthHist;
my %typeHist;
my %repeats;
my %repeatCat;
my %repeatClass;
my %multiNoRepeat;
my %multiNoRepeatTyped;
my %ids;
my %unique;
my %noHit;
my %proteins;
my %snos;
my %miRNAs;
my %lincs;
my %antisense;
my %catchall;
my %mito;
my %ribo;
my %tRNA;
my %mtRibo;
my %mttRNA;
my %uniqueRep;
my %ambig;
my %strands;
my %chrom;
my %ambigMap;
my %misMatch;
my %cntHist;

my $gfh = FileUtil::openFileHandle($gtf);
while (<$gfh>) {
#last;
	chomp;
	my @d = split /\t/;
	next unless $d[2] eq "gene";
	$d[8] =~ /gene_id .(\S+)\S\S/;
	$strands{$1} = $d[6];
	$chrom{$1} = $d[0];
}
$gfh->close();

my @order;
my @altOrder;
my $ofh = FileUtil::openFileHandle($order);
while (<$ofh>) {
	chomp;
	my @d = split /\t/;
	if ($d[1]) {
		push @order,$d[0];
		push @altOrder,$d[1];
	} else {
		push @order,$d[0];
		push @altOrder,$d[0];
	}
}
$ofh->close();

my $lfh = FileUtil::openFileHandle($list);
while (<$lfh>) {
	chomp;
	my $f = $_;
	print STDERR $f,"\n";
#	print $f,"\n";
	my $id;
#	if ($f =~ /(test.\w+.data)$/) {
	if ($f =~ /(\w+).sam.gz$/) {
		$id = $1;
	} else {
		print STDERR "ERROR: Could not identify id in file $f\n";
		next;
	}
#	print $id,"\n";
	$ids{$id} = 1;
	my $fh = FileUtil::openFileHandle("$f");
#	$xxx = 0;
	while (<$fh>) {
#print $_;
#		last if $xxx > 10;
		chomp;
#		next unless /unique/;
                next if($_ =~ /^\@/o);
		my $l = $_;
		$totalReads{$id}++;
		my $type = 0;
		if ($l =~ /XT:i:(\d+)/) {
			$type = $1;
		} else {
			print STDERR "ERROR: No type found $_\n";
			next;
		}
		my @d = split /\t/;
		my $len = length($d[9]); 
		$lengthHist{Total}{$id}{$len}{$type}++;
		$typeHist{$id}{$type}++;
		my $multi = 0;
		if ($l =~ /multi_match/) {
			$multi = 1;
			$lengthHist{Multi}{$id}{$len}{$type}++;
			$multiNoRepeat{$id}++;
			$multiNoRepeatTyped{$id}{$type}++;
		} else {
			$lengthHist{Single}{$id}{$len}{$type}++;
		}
		if ($d[1] == 4) {
			### no match found
			$noHit{$id}++;
			next;
		}
		if ($type >  1) {
			next;
		}
		if ($l =~ /mis_match/) {
			if ($l =~ /XL:Z:/) {
			} else {
				$misMatch{$id}++;
				next;
			}
		}
		if ($d[2] eq "chrM") {
			my $moveOn = 0;
			if ($l =~ /NH:i:(\d+)/) {
				if ($1 == 1) {
					$unique{$id}{mito}++;
					$moveOn = 1;
				}
			}
			if ($l =~ /Mt_rRNA/) {
				if ($l =~ /ambiguous/) {
					if (testSimilarity($l,"==","rRNA")) {
						$mtRibo{$id}++;
						next;
					}
				} else {
					$mtRibo{$id}++;
					next;
				}
			}
			if ($l =~ /Mt_tRNA/) {
				if ($l =~ /ambiguous/) {
					if (testSimilarity($l,"==","tRNA")) {
						$mttRNA{$id}++;
						next;
					}
				} else {
					$mttRNA{$id}++;
					next;
				}
			}
			if ($l =~ /XG:Z:.*unique_match.(\S+?);/) {
				$mito{$id}{$1}++;
				$moveOn = 1;
			}
			next if $moveOn;
		}
		if ($l =~ /XR:Z:(\S+)\tXL:Z:(\S+)/) {
			my $rep = $1;
			my $rcl = $2;
			my $good = 0;
			my $v = "~~";
			if ($l =~ /NH:i:1\s/ && $l !~ /mis_match/ && $l !~ /XL:Z:tRNA/) {
				if ($l =~ /XG:Z:(\S+)/) {
					my $tmp = $1;
					my @set = split /==/;
					foreach my $i (@set) {
						my @set2 = split /$v/,$i;
						$good++ if $set2[1];
					}
				}
			}
			if ($good) {
			
			} else {
				if ($d[2] eq "chrR") {
					$rep = "Ribo";
					$rcl = "Ribo";
				}
				unless ($type == 2) {
					### is a repeat
					$repeats{$id}++;
					$repeatCat{$id}{$rep}++;
					$repeatClass{$id}{$rcl}++;
					$repeatCat{all}{$rep}++;
					$repeatClass{all}{$rcl}++;
					$l =~ /NH:i:(\d+)/;
					my $cnt = $1;
					$cnt = 10 if $cnt > 10;
					$uniqueRep{$id}[$cnt]++;
					if ($rcl eq "rRNA" or $rcl eq "Ribo") {
						$ribo{$id}++;
	#					print $l,"\t$rep\t$rcl\n";
						next;
					}
					if ($rcl eq "tRNA") {
						$tRNA{$id}++;
						next;
					}
					next if $cnt > 1;
					next if $l =~ /mis_match/;
					next if $d[1] == 4;
				}
			}
		}
		if ($multi) {
		} elsif ($type == 1) {
			if ($l =~ /ambiguous_match.(\S+)/) {
				### ambiguous (matches 2 or more genes)
#				print "XXX\t$l\n";
				my $hits = $1;
				$ambig{$id}++;
				my @hits = split /==/,$hits;
				my @goodHits;
				my %goodIds;
				foreach my $i (@hits) {
#					print "VVV $i\n";
					my @d = split /~~/,$i;
					if ($d[1]) {
						if ($d[0] =~ /^(\S+?);/) {
							push @goodHits,$i unless exists($goodIds{$1});
#							my $v = @goodHits;
#							print "WWW $v $1 $d[1]\n";
							$goodIds{$1} = 1;
						}
					}
				}
				my $goodCnt = @goodHits;
				if ($goodCnt) {
					my %ambigHits;
					my %seen;
					foreach my $i (@goodHits) {
#						print "YYY $i\n";
						if ($i =~ /(\S+);.*protein_coding/) {
							### gene
							my $pid = $1;
							$proteins{$id}{$pid}++;
							$unique{$id}{protein}++ unless exists($seen{protein});
							$seen{protein} = 1;
							$ambigHits{$pid}++ if $goodCnt > 1;
#							$xxx++;
#							print "ooo\t$id\t$pid\t$i\n" unless exists($chrom{$pid});
						}
						elsif ($i =~ /(\S+);.*snoRNA/) {
							### snoRNA
							my $pid = $1;
	#						print "$id\t$pid\t$i\n";
							$snos{$id}{$pid}++;
							$unique{$id}{snoRNA}++ unless exists($seen{snoRNA});
							$seen{snoRNA} = 1;
							$ambigHits{$pid}++ if $goodCnt > 1;
						}
						elsif ($i =~ /(\S+);.*miRNA/) {
							### miRNA
							my $pid = $1;
							$miRNAs{$id}{$pid}++;
							$unique{$id}{miRNA}++ unless exists($seen{miRNA});
							$seen{miRNA} = 1;
							$ambigHits{$pid}++ if $goodCnt > 1;
#							print "$id\t$pid\t$i\t$unique{$id}{miRNA}\t$finished\n";
						}
						elsif ($i =~ /(\S+);.*lincRNA/) {
							### lincRNA
							my $pid = $1;
	#						print "$id\t$pid\t$i\n";
							$lincs{$id}{$pid}++;
							$unique{$id}{lincRNA}++ unless exists($seen{lincRNA});
							$seen{lincRNA} = 1;
							$ambigHits{$pid}++ if $goodCnt > 1;
						}
                                                elsif ($i =~ /(\S+);.*antisense/) {
                                                        ### antisense
                                                        my $pid = $1;
        #                                               print "$id\t$pid\t$i\n";
                                                        $antisense{$id}{$pid}++;
                                                        $unique{$id}{antisense}++ unless exists($seen{antisense});
                                                        $seen{antisense} = 1;
                                                        $ambigHits{$pid}++ if $goodCnt > 1;
                                                }
                                                elsif ($i =~ /(\S+);\d+\.\d+\.(\+|\-)\~\~(\S+)\~\~/) {
                                                        ### catchall
                                                        # ENSG00000242590.1;990413.990595.+~~sense_intronic~~NOVEL~~RP11-54O7.14
                                                        my $pid = $1;
        #                                               print "$id\t$pid\t$i\n";
                                                        $catchall{$id}{$pid}++;
                                                        $unique{$id}{$3}++ unless exists($seen{$3});
                                                        $seen{$3} = 1;
                                                        $ambigHits{$pid}++ if $goodCnt > 1;
                                                }
						else {
							$i =~ /^(\S+);/;
							my $pid = $1;
							$ambigHits{$pid}++ if $goodCnt > 1;
						}
					}
					foreach my $i (keys %ambigHits) {
						foreach my $j (keys %ambigHits) {
							next if $i eq $j;
#							print "ZZZ\t$id\t$i\t$j\n";#ENSG00000185989.9
							$ambigMap{$id}{$i}{$j} = 1;
						}
					}
				}
			}
			elsif ($d[2] eq "chrM") {
				### mitochondrial
				if ($l =~/XG:Z:.*unique_match.(\S+?);/) {
					$mito{$id}{$1}++;
				}
				$unique{$id}{mito}++;
			}
			elsif ($l =~ /XG:Z:\S+\!(\S+?);.*protein_coding/) {
				### gene
#				$xxx++;
				my $pid = $1;
				$proteins{$id}{$pid}++;
				$unique{$id}{protein}++;
#				print "RRR $id\t$pid\t$l\n" unless exists($chrom{$pid});
			}
			elsif ($l =~ /XG:Z:\S+\!(\S+?);.*snoRNA/) {
				### snoRNA
				my $pid = $1;
				$snos{$id}{$pid}++;
				$unique{$id}{snoRNA}++;
			}
			elsif ($l =~ /XG:Z:\S+\!(\S+?);.*miRNA/) {
				### miRNA
				my $pid = $1;
				$miRNAs{$id}{$pid}++;
				$unique{$id}{miRNA}++;
			}
			elsif ($l =~ /XG:Z:\S+\!(\S+?);.*lincRNA/) {
				### lincRNA
				my $pid = $1;
				$lincs{$id}{$pid}++;
				$unique{$id}{lincRNA}++;
			}
                        elsif ($l =~ /XG:Z:\S+\!(\S+?);.*antisense/) {
                                ### antisense
                                my $pid = $1;
                                $antisense{$id}{$pid}++;
                                $unique{$id}{antisense}++;
                        }
                        elsif ($l =~ /XG:Z:\S+\!(\S+?);.*(polymorphic_pseudogene|processed_transcript|pseudogene|sense_intronic|sense_overlapping|misc_RNA)/) {
                                # GENE:seq303312!chr1!chr1!Exon!unique_match!ENSG00000188157.9;990204.991496.+~~protein_coding~~KNOWN~~AGRN==ENSG00000188157.9;990204.991492.+~~protein_coding~~KNOWN~~AGRN
                                ### catchall
                                my $pid = $1;
                                $catchall{$id}{$pid}++;
                                $unique{$id}{misc}++;
                        }
			elsif ($l =~ /XG:Z:\S+no_match/) {
				### reads that no known gene and are not repetitive
				$unique{$id}{nomatch}++;
			} else {
				### other pseudogene, other RNAs, mismatches
				$unique{$id}{misc}++;
			}
		}
	}
	$fh->close();
}
$lfh->close();

my $sfho = new FileHandle(">$summaryOut");
my $l = join "\t","",@altOrder;
print $sfho $l,"\n";
### total hits;
my @l;
push @l,"All Reads";
foreach my $i (@order) {
	push @l,0+$totalReads{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### no hits;
my @l;
push @l,"No Hits";
foreach my $i (@order) {
	push @l,0+$noHit{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 0;
my @l;
push @l,"Type 0";
foreach my $i (@order) {
	push @l,0+$typeHist{$i}{0};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 1;
my @l;
push @l,"Type 1";
foreach my $i (@order) {
	push @l,0+$typeHist{$i}{1};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 2;
my @l;
push @l,"Type 2";
foreach my $i (@order) {
	push @l,0+$typeHist{$i}{2};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 3;
my @l;
push @l,"Type 3";
foreach my $i (@order) {
	push @l,0+$typeHist{$i}{3};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 1;
my @l;
push @l,"Multi Type 1";
foreach my $i (@order) {
	push @l,0+$multiNoRepeatTyped{$i}{1};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 2;
my @l;
push @l,"Multi Type 2";
foreach my $i (@order) {
	push @l,0+$multiNoRepeatTyped{$i}{2};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### type 3;
my @l;
push @l,"Multi Type 3";
foreach my $i (@order) {
	push @l,0+$multiNoRepeatTyped{$i}{3};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### ambiguous matches
my @l;
push @l,"Ambiguous";
foreach my $i (@order) {
	push @l,0+$ambig{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### mis matches
my @l;
push @l,"MisMatches";
foreach my $i (@order) {
	push @l,0+$misMatch{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### unique repeats
for (my $i = 1; $i <= 10; $i++) {
	my @l;
	push @l,"Unique Rep $i";
	foreach my $j (@order) {
		push @l,0+$uniqueRep{$j}[$i];
	}
	my $l = join "\t",@l;
	print $sfho $l,"\n";
}

### repetitive;
my @l;
push @l,"Repetitive";
foreach my $i (@order) {
	push @l,0+$repeats{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### ribosomal;
my @l;
push @l,"Ribosomal";
foreach my $i (@order) {
	push @l,0+$ribo{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### tRNA;
my @l;
push @l,"tRNA";
foreach my $i (@order) {
	push @l,0+$tRNA{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### multi hits;
my @l;
push @l,"Multiple Hits";
foreach my $i (@order) {
	push @l,0+$multiNoRepeat{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### no match;
my @l;
push @l,"No Match";
foreach my $i (@order) {
	push @l,0+$unique{$i}{nomatch};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### misc
my @l;
push @l,"Misc";
foreach my $i (@order) {
	push @l,0+$unique{$i}{misc};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### Mito;
my @l;
push @l,"Mitochondria";
foreach my $i (@order) {
	push @l,0+$unique{$i}{mito};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### mito ribosomal;
my @l;
push @l,"Mito Ribosomal";
foreach my $i (@order) {
	push @l,0+$mtRibo{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### mttRNA;
my @l;
push @l,"mito tRNA";
foreach my $i (@order) {
	push @l,0+$mttRNA{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### mitos;
my @l;
push @l,"Number of Mito Proteins";
foreach my $i (@order) {
	my $cnt;
	foreach my $j (keys %{$mito{$i}}) {
		$cnt++ if $mito{$i}{$j} > 1;
	}
	push @l,$cnt;
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### protein coding;
my @l;
push @l,"Protein Coding";
foreach my $i (@order) {
	push @l,0+$unique{$i}{protein};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### proteins;
my %plus;
my %minus;
my %chrPlace;
my $chroms;
my @l;
push @l,"Number of Proteins";
foreach my $i (@order) {
	my $cnt;
	foreach my $j (keys %{$proteins{$i}}) {
		push @{$cntHist{protein}{$i}},$proteins{$i}{$j};
		if ($proteins{$i}{$j} > 1) {
			$cnt++;
			if ($strands{$j} eq "+") {
				$plus{$i}++;
			} else {
				$minus{$i}++;
			}
			$chrPlace{$i}{$chrom{$j}}++;
			$chroms{$chrom{$j}}++;
		}
	}
	push @l,$cnt;
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### protein plus strand;
my @l;
push @l,"Plus Strand";
foreach my $i (@order) {
	push @l,0+$plus{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### protein minus strand;
my @l;
push @l,"Minus Strand";
foreach my $i (@order) {
	push @l,0+$minus{$i};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### snoRNA
my @l;
push @l,"snoRNA";
foreach my $i (@order) {
	push @l,0+$unique{$i}{snoRNA};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### snoRNAs;
my @l;
push @l,"Number of snoRNAs";
foreach my $i (@order) {
	my $cnt;
	foreach my $j (keys %{$snos{$i}}) {
		push @{$cntHist{snoRNA}{$i}},$snos{$i}{$j};
		$cnt++ if $snos{$i}{$j} > 1;
	}
	push @l,$cnt;
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### miRNA
my @l;
push @l,"miRNA";
foreach my $i (@order) {
	push @l,0+$unique{$i}{miRNA};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### miRNAs;
my @l;
push @l,"Number of miRNAs";
foreach my $i (@order) {
	my $cnt;
	foreach my $j (keys %{$miRNAs{$i}}) {
		$cnt++ if $miRNAs{$i}{$j} > 1;
		push @{$cntHist{miRNA}{$i}},$miRNAs{$i}{$j};
	}
	push @l,$cnt;
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### lincRNA
my @l;
push @l,"lincRNA";
foreach my $i (@order) {
	push @l,0+$unique{$i}{lincRNA};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### lincs;
my @l;
push @l,"Number of lincRNAs";
foreach my $i (@order) {
	my $cnt;
#	my $cnt = keys %{$lincs{$i}};
	foreach my $j (keys %{$lincs{$i}}) {
		$cnt++ if $lincs{$i}{$j} > 1;
		push @{$cntHist{linc}{$i}},$lincs{$i}{$j};
	}
	push @l,$cnt;
}
my $l = join "\t",@l;
print $sfho $l,"\n";
### antisense
my @l;
push @l,"antisense";
foreach my $i (@order) {
        push @l,0+$unique{$i}{antisense};
}
my $l = join "\t",@l;
print $sfho $l,"\n";

### antisense;
my @l;
push @l,"Number of antisense";
foreach my $i (@order) {
        my $cnt;
#       my $cnt = keys %{$antisense{$i}};
        foreach my $j (keys %{$antisense{$i}}) {
                $cnt++ if $antisense{$i}{$j} > 1;
                push @{$cntHist{antisense}{$i}},$antisense{$i}{$j};
        }
        push @l,$cnt;
}

my $l = join "\t",@l;
print $sfho $l,"\n";



print $sfho "\n";
my $l = join "\t","",@altOrder;
print $sfho $l,"\n";
my $cnt = 0;
foreach my $j (sort { $repeatCat{all}{$b} <=> $repeatCat{all}{$a};} keys %{$repeatCat{all}}) {
	### repeat categories
	my @l;
	push @l,"RCat:$j";
	foreach my $i (@order) {
		push @l,0+$repeatCat{$i}{$j};
	}
	my $l = join "\t",@l;
	print $sfho $l,"\n";
	$cnt++;
	last if $cnt == 10;
}

print $sfho "\n";
my $l = join "\t","",@altOrder;
print $sfho $l,"\n";
my $cnt = 0;
foreach my $j (sort { $repeatClass{all}{$b} <=> $repeatClass{all}{$a};} keys %{$repeatClass{all}}) {
	### repeat categories
	my @l;
	push @l,"RClass:$j";
	foreach my $i (@order) {
		push @l,0+$repeatClass{$i}{$j};
	}
	my $l = join "\t",@l;
	print $sfho $l,"\n";
	$cnt++;
	last if $cnt == 10;
}

print $sfho "\n";
my $l = join "\t","",@altOrder;
print $sfho $l,"\n";
foreach my $j (sort keys %chroms) {
	### genes by chrom
	my @l;
	push @l,$j;
	foreach my $i (@order) {
		push @l,0+$chrPlace{$i}{$j};
	}
	my $l = join "\t",@l;
	print $sfho $l,"\n";
}


$sfho->close();

my $lfho = new FileHandle(">$baseHistOut.length.vs.type.hist");
#		$lengthHist{Total}{$id}{$len}{$type}++;
my @l;
push @l,"","No Match","Human","Bacterial","Mixed";
my $l = join "\t",@l;
print $lfho $l,"\n";
foreach my $i (@order) {
	for (my $j = 17; $j <= 51; $j++) {
		my @l;
		push @l,$i,$j;
		for (my $k = 0; $k < 4; $k++) {
			push @l,(0+$lengthHist{Total}{$i}{$j}{$k});
		}
		my $l = join "\t",@l;
		print $lfho $l,"\n";
	}
}
$lfho->close();

	#		$lengthHist{Total}{$id}{$len}{$type}++;
foreach my $i (keys %ambigMap) {
	my $lfho = new FileHandle(">$baseHistOut.$i.ambiguous.links");
#	print $lfho "AAA $i\n";
	foreach my $j (keys %{$ambigMap{$i}}) {
#		print $lfho "BBB $j\n";
		foreach my $k (keys %{$ambigMap{$i}{$j}}) {
#			print $lfho "CCC $k\n";
			print $lfho "$j\t$k\t1\n";
		}
	}
	$lfho->close();
}

my $lfho = new FileHandle(">$baseHistOut.abundance.cnts");
foreach my $i (keys %cntHist) {
	foreach my $j (keys %{$cntHist{$i}}) {
		my $l = join ",",@{$cntHist{$i}{$j}};
		print $lfho "$i\t$j\t$l\n";
	}
}
$lfho->close();

my $lfho = new FileHandle(">$baseHistOut.gene.cnts");
{
	my %tmp;
	foreach my $i (keys %proteins) {
		foreach my $j (keys %{$proteins{$i}}) {
			$tmp{$j}{$i} = $proteins{$i}{$j};
		}
	}
	foreach my $i (keys %miRNAs) {
		foreach my $j (keys %{$miRNAs{$i}}) {
			$tmp{$j}{$i} = $miRNAs{$i}{$j};
		}
	}
	foreach my $i (keys %lincs) {
		foreach my $j (keys %{$lincs{$i}}) {
			$tmp{$j}{$i} = $lincs{$i}{$j};
		}
	}
        foreach my $i (keys %antisense) {
                foreach my $j (keys %{$antisense{$i}}) {
                        $tmp{$j}{$i} = $antisense{$i}{$j};
                }
        }
	foreach my $i (keys %snos) {
		foreach my $j (keys %{$snos{$i}}) {
			$tmp{$j}{$i} = $snos{$i}{$j};
		}
	}
        foreach my $i (keys %catchall) {
                foreach my $j (keys %{$catchall{$i}}) {
                        $tmp{$j}{$i} = $catchall{$i}{$j};
                }
        }
	my $l = join "\t","",@altOrder;
	print $lfho "$l\n";
	foreach my $i (keys %tmp) {
		my @set;
		push @set,$i;
		foreach my $j (@order) {
			push @set,0+$tmp{$i}{$j};
		}
		my $l = join "\t",@set;
		print $lfho "$l\n";
	}
}
$lfho->close();

sub testSimilarity {
	my ($str,$sep,$query) = @_;
	
	my @pieces = split /$sep/,$str;
	my $good = 1;
	foreach my $i (@pieces) {
		if ($i =~ /$query/) {
		} else {
			$good = 0;
		}
	}
	
	return $good;
}



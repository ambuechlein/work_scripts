#!/usr/bin/env perl
use strict;
use warnings;

my $gff = shift;
my $chr = shift;
my $bamList = shift;
my $dir = shift;
my $mod = shift;
my $error = shift;
my $min = shift;
$error = 2 unless defined($error);
$min = 10 unless defined($min);
$dir = "." unless defined($dir);
my $first = 0;
warn $chr,"\n";
open(my $gfh, $gff) or die "Can't open $gff: $!\n";
my %genes;
my %geneName;
my %geneChr;
my %geneStatus;
my %geneType;
my %exons;
warn "Parsing GFF: $gff\n";
while (<$gfh>) {
  my $l = $_;
  chomp;
  my @d = split /\t/;
  next unless $d[0] eq $chr;
  if ($d[2] eq "gene") {
    my $b = $d[3];
    my $e = $d[4];
    my $str = $d[6];
    $str = $str eq "+" ? 0:1;
    my $geneId;
    if ($d[8] =~ /gene_id\s+\"(\S+?)\"/ || $d[8] =~ /ID=(.+?);/) {
      $geneId = $1;
    } else {
      warn "ERROR: no gene id found $l\n";
      exit;
    }
    if ($d[8] =~ /gene_type\s+\"(\S+?)\"/ || $d[8] =~ /Name=(.+?);/) {
      $geneType{$geneId} = $1;
    } else {
      warn "ERROR: no gene type found $l\n";
    }
    if ($d[8] =~ /gene_status\s+\"(\S+?)\"/ || $d[8] =~ /locus_tag=(\S+)/ || $d[8] =~ /Name=(.+?);/) {
      $geneStatus{$geneId} = $1;
    } else {
      if($d[8] =~ /(level\s*\d+)\;/){
         $geneStatus{$geneId} = $1;
      }else{
        warn "ERROR: no gene status found $l\n";
      }
    }
    if ($d[8] =~ /gene_name\s+\"(\S+?)\"/ || $d[8] =~ /Name=(.+?);/) {
      $geneName{$geneId} = $1;
    } else {
      warn "ERROR: no gene name found $l\n";
    }
    $geneChr{$geneId} = $d[0];
    for (my $i = $b; $i <= $e; $i++) {
      $genes{$d[0]}{$str}{$i}{$geneId} = 1
    }
  } elsif ($d[8] =~ /gene_id/) {
    if ($d[2] eq "exon") {
      my $b = $d[3];
      my $e = $d[4];
      my $str = $d[6] eq "+" ? 0:1;
      my $geneId;
      if ($d[8] =~ /gene_id\s+\"(\S+?)\"/) {
        $geneId = $1;
        if($d[1] eq 'hg19_rmsk' or $d[1] eq 'hg38_rmsk'){
          $geneType{$geneId} = 'repeat';
          $geneStatus{$geneId} = 'RepeatMasker';
          $geneName{$geneId} = $geneId;
        }
      } else {
        warn "ERROR: no gene id found $l\n";
      }
      my $exonId = "$geneId;$b.$e.$d[6]";
      for (my $i = $b; $i <= $e; $i++) {
        $exons{$d[0]}{$str}{$i}{$exonId} = 1;
      }
    }
  } else {
    if ($d[2] eq "CDS" || $d[2] eq "exon") {
      my $b = $d[3];
      my $e = $d[4];
      my $str = $d[6] eq "+" ? 0:1;
      my $geneId;
      if ($d[8] =~ /Parent=(\S+?);/) {
        $geneId = $1;
      } else {
        warn "ERROR: no gene id found $l\n";
      }
      if(not defined $geneType{$geneId}){
#      $genes{$d[0]}{$str}{$i}{$geneId} = 1
        my $p;
        my @kb = keys %{$genes{$d[0]}{$str}{$b}};
        my @ke = keys %{$genes{$d[0]}{$str}{$e}};
        my $c;
        if(scalar @kb < 1){
#ID=id38;Parent=rna37;gbkey=tRNA;product=tRNA-Arg
          my @test = split(/\;/, $d[8]);
#          warn join("\n", @test), "\n\n";
          foreach my $att (@test){
            my ($n, $v) = split(/\=/, $att);
            $geneType{$geneId} = $v if($n eq 'gbkey');
            $geneStatus{$geneId} = $v if($n eq 'gbkey');
            $geneName{$geneId} = $v if($n eq 'product');
          }
        }
        if(scalar @kb >= 1 and scalar @ke >= 1){
          foreach my $a (@kb){
            foreach my $b (@ke){
              if($a eq $b){
                $p = $a;
                $c++;
              }
            }
          } 
          warn "Warning: Too may gene keys to set status for $geneId: ",join(" : ", @kb), " | ", join(" : ", @ke),  "\n" if($c > 1);
          warn "Randomly selected $p for annotation purposes\n" if($c > 1); 
#          my @kb = keys %{$genes{$d[0]}{$str}{$b}};
#          my $p = shift @k;
          $geneType{$geneId} = $geneType{$p};
          $geneStatus{$geneId} = $geneStatus{$p};
          $geneName{$geneId} = $geneName{$p};
        }
      }
      my $exonId = "$geneId;$b.$e.$d[6]";
      for (my $i = $b; $i <= $e; $i++) {
        $exons{$d[0]}{$str}{$i}{$exonId} = 1;
      }
    }
  }
}
close $gfh;

warn "Stage 1 $chr done\n";
my %seen_base;
open(my $blfh, $bamList) or die "Can't open $bamList: $!\n";
while (<$blfh>) {
  chomp;
  my $bam = $_;
  warn "$bam\n";
  my $bfh;
  my $cntr = 0;
  while (1) {
    if ($bam =~ /\.sam$/ || $bam =~ /\.sam.gz/) {
       open($bfh, $bam =~ /.gz(ip)?$/ ? "zcat $bam |" : $bam =~ /.bz(ip)?2$/ ? "bzcat $bam |" : $bam) || die("Open error: $bam");
    } else {
       open($bfh, "samtools view $bam |");
    }
    last if defined($bfh);
    sleep 5;
    $cntr++;
    if ($cntr == 20) {
      warn "ERROR: Could not open file $bam\n";
      die "ERROR: Could not open file $bam\n";
      exit 1;
    }
  }
  my $baseName;
#  if ($bam =~ /(\w+)\.fastq\.nonredundant\.fasta_thout.*sorted.[s|b]am/) {
  if ($bam =~ /((\w|-)+)\.fastq(\.gz)?\.nonredundant\.fasta_(thout|hisat).*sorted.[s|b]am/) {
    $baseName = $1;
    warn "------------------------ WARNING ------------------------\n$baseName used twice\n------------------------ WARNING ------------------------\n" if(defined $seen_base{$baseName}); 
    $seen_base{$baseName} = 1;
  } else {
    warn "ERROR: unable to parse file name $bam\n";
    die "ERROR: unable to parse file name $bam\n";
    exit 1;
  }
  open(my $fo, ">$dir/$baseName.$chr.$mod.counts") or die "Can't open $dir/$baseName.$chr.$mod.counts for writing: $!\n";
  while (<$bfh>) {
    my $l = $_;
    next if /^\@/;
    chomp;
    my @d = split /\t/;
    next if $d[1] == 4;
    next if $d[2] ne $chr; # Aaron moved up
    my $str = $d[1] & 16 ? 1 : 0;
    if (/NM:i:(\d+)/ && $1 > $error) {
      print $fo "$d[0]\t$d[2]\t$chr\tGene\terroneous_match\n";
      next;
    }
    if (/NH:i:(\d+)/ && $1 > 1) {
      print $fo "$d[0]\t$d[2]\t*\tGene\tmulti_match\n";
      next;
    }
#    next if $d[2] ne $chr; # Aaron commented out
    my @cigar;
    while ($d[5] =~ s/(\d+)(\D)//) {
      my $c = $1;
      my $v = $2;
      my $r = $v x $c;
      push @cigar,(split //,$r);
    }
    my $p = $d[3];
    my $len = 0;
    my %geneList;
    my %geneList2;
    my %exonList;
    while (@cigar) {
      my $action = shift @cigar;
      my @genes = ();
      my @exons = ();
      if ($action eq "M") {
        if(defined $genes{$d[2]}{$str}{$p}){
          @genes = keys %{$genes{$d[2]}{$str}{$p}};
        }

        foreach my $j (@genes) {
          $geneList{$j}++;
        }
        if(defined $exons{$d[2]}{$str}{$p}){
          @exons =  keys %{$exons{$d[2]}{$str}{$p}};
        }
        my %first;
        foreach my $j (@exons) {
          $exonList{$j}++;
          $j =~ /(\S+);/;
          $geneList2{$1}++ unless exists($first{$1});
          $first{$1}++;
        }
        $p++;
        $len++;
      } elsif ($action eq "N") {
        $p++;
      } elsif ($action eq "I") {
        $p++;
        $len++;
      } elsif ($action eq "D") {
        $len++;
      }
    }
    my $numGenes = keys %geneList;
    if ($numGenes > 0) {
      my $good = 0;
      my $pick;
      my @list;
      foreach my $i (sort { $geneList{$b} <=> $geneList{$a};} keys %geneList) {
        if ($good) {
          if ($geneList{$i} + $error >= $len) {
            $good++;
            push @list,$i;
          }
        } else {
          $good++ if $geneList{$i} + $error >= $len;
          push @list,$i;
          $pick = $i;
        }
      }
      if ($good == 1) {
        print $fo "$d[0]\t$d[2]\t$chr\tGene\tunique_match\t$pick\t$geneType{$pick}\t$geneStatus{$pick}\t$geneName{$pick}\n";
      } elsif ($good > 1) {
        my @set;
        foreach my $i (@list) {
          push @set,join "~~",$i,$geneType{$i},$geneStatus{$i},$geneName{$i};
        }
        my $list = join "==",@set;
        print $fo "$d[0]\t$d[2]\t$chr\tGene\tambiguous_match\t$list\n";
      } else {
        my $list = join ",",@list;
        print $fo "$d[0]\t$d[2]\t$chr\tGene\tmis_match\t$list\n";
      }
    } else {
      print $fo "$d[0]\t$d[2]\t$chr\tGene\tno_match\n";
    }
    my $numExons = keys %geneList2;
    if ($numExons > 0) {
      my $good = 0;
      my $pick;
      my %list;
      foreach my $i (sort { $geneList2{$b} <=> $geneList2{$a};} keys %geneList2) {
        if ($good) {
          if ($geneList2{$i} + $error >= $len) {
            $good++;
            $list{$i}++;
          }
        } else {
          $good++ if $geneList2{$i} + $error >= $len;
          $pick = $i;
          $list{$i}++;
        }
      }
      my @list;
      foreach my $i (keys %exonList) {
        $i =~ /(\S+);/;
        my $gId = $1;
        push @list,$i if (exists($list{$gId}));
      }
      my @set;
      foreach my $i (@list) {
        $i =~ /(\S+);/;
if(not defined $geneName{$1} or not defined $geneType{$1} or not defined $geneStatus{$1}){
 warn "-----------Aaron-----------\n";
 warn "key: $1\ni: $i\ngeneType: $geneType{$1}\ngeneStatus: $geneStatus{$1}\ngeneName: $geneName{$1}\n";
 warn "---------------------------\n\n";
 die if(not defined $geneName{$1} or not defined $geneType{$1} or not defined $geneStatus{$1});
}
        push @set,join "~~",$i,$geneType{$1},$geneStatus{$1},$geneName{$1};
      }
      my $list = join "==",@set;
      if ($good == 1) {
        $pick =~ /(\S+);/;
        print $fo "$d[0]\t$d[2]\t$chr\tExon\tunique_match\t$list\n";
      } elsif ($good > 1) {
        print $fo "$d[0]\t$d[2]\t$chr\tExon\tambiguous_match\t$list\n";
      } else {
        print $fo "$d[0]\t$d[2]\t$chr\tExon\tmis_match\t$list\n";
      }
    } else {
      print $fo "$d[0]\t$d[2]\t$chr\tExon\tno_match\n";
    }
  }
  close $bfh;
  warn "Fini\n";
  close $fo;
}
close $blfh;
$first++;


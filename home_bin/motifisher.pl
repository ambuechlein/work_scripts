#!/usr/bin/env perl -w
use strict;
use File::Basename;
use Getopt::Std;
use POSIX;  # ceil()
my %opts;
getopts('s:p:t:l', \%opts);
@ARGV == 3 or die "Usage: " . basename($0) . " [-lpst] <mats> <myseqs.fasta> <BGseqs.fasta>\n";
my $pthresh = exists $opts{'t'} ? $opts{'t'} : 1;
my $bg_seq_percent = exists $opts{'s'} ? $opts{'s'}/100 : 0.1;
# Default value of 10% is roughly the criterion of Elkon et al.
my ($matfile, $fgfile, $bgfile) = @ARGV;
my $clover = "./clover -u0";
$clover .= " -p $opts{'p'}" if exists $opts{'p'};
$clover .= " -l" if exists $opts{'l'};

my @bgseqlens = get_seqlens($bgfile);
my @fgseqlens = get_seqlens($fgfile);

my %matlens;
open FILE, $matfile or die $!;
my $matname;
my $matlen = 0;
while (<FILE>) {
    if ( /^\s*>(.*?)\s*$/ ) {
	$matlens{$matname} = $matlen if defined $matname;
	$matname = $1;
	$matlen = 0;
    } else {
	++$matlen;
    }
}
$matlens{$matname} = $matlen if defined $matname;  # don't forget the last one!

my %bgseqscores;
my %topscores;

my @bgoutput=split/\n+/,`$clover $matfile $bgfile 2>/dev/null`;
for (@bgoutput) {
    if ( /^\s*>/ ) {
	push @{$bgseqscores{$_}}, $topscores{$_} for keys %topscores;
	%topscores = ();
    } elsif ( /^(.*?)\s+\d+\s+-\s+\d+\s+[+-]\s+\w+\s+(\S+)/ ) {
	$topscores{$1} = $2
	    unless exists $topscores{$1} and $topscores{$1} > $2;
    }
}
push @{$bgseqscores{$_}}, $topscores{$_} for keys %topscores;  # !

my $bgseqnum = @bgseqlens;

#0.1 approx criterion from Elkon et al.
#used floor() for the results in the Clover NAR 2004 paper:
my $tsn = POSIX::ceil($bgseqnum * $bg_seq_percent);

my %seqthresh;
my %p;
for (keys %matlens) {
    if(!$bgseqscores{$_}){next}
    $seqthresh{$_} = (sort { $b <=> $a } @{$bgseqscores{$_}})[$tsn - 1];
}

my %fgseqs;
my %thisseq;
my @fgoutput=split/\n+/,`$clover $matfile $fgfile 2>/dev/null`;
for (@fgoutput) {
    if ( /^\s*>/ ) {
	++$fgseqs{$_} for keys %thisseq;
	%thisseq = ();
    } elsif ( /^(.*?)\s+\d+\s+-\s+\d+\s+[+-]\s+\w+\s+(\S+)/ ) {
	$thisseq{$1} = 1 if $2 >= $seqthresh{$1};
    }
}
++$fgseqs{$_} for keys %thisseq;  # !

my %bgseqs;
%thisseq = ();
for (@bgoutput) {
    if ( /^\s*>/ ) {
	++$bgseqs{$_} for keys %thisseq;
	%thisseq = ();
    } elsif ( /^(.*?)\s+\d+\s+-\s+\d+\s+[+-]\s+\w+\s+(\S+)/ ) {
	$thisseq{$1} = 1 if $2 >= $seqthresh{$1};
    }
}
++$bgseqs{$_} for keys %thisseq;  # !

my @lookup = (0);  # needed by fisher
for (keys %matlens) {
    $fgseqs{$_} = 0 unless defined $fgseqs{$_};
$p{$_}=sprintf"%6.4g",fisher($fgseqs{$_}, @fgseqlens - $fgseqs{$_},$bgseqs{$_}, @bgseqlens - $bgseqs{$_});
}

print "Motifish: Motif overrepresentation by Fisher Exact tests\n";
my $overflag;my %output;
for(@fgoutput)
 {if ( /Sorry, couldn.t understand the matrix file/ ) {print "Motifisher failed: $_";exit}
 elsif(/Clover: Cis-eLement OVERrepresentation|^Compiled on/){next}
 elsif(/^(Motif +)Raw score.*/){print "${1}P-value from Fisher Exact Tests\n";$overflag=1}
 elsif(/^(\*+ +Motif Instances with Score)/i){print "$1 above Threshold\n";$overflag=0;}
 elsif($overflag&&/^(.*?\S)(\s+)[\-\.\d]+\s+$/){if($p{$1}<=$pthresh){print "$1$2$p{$1}\n";$output{$1}=1;}}
 elsif(/^\s*(.*?)\s+(\d+)\s+\-\s+(\d+)\s+([+-])\s+(\w+)\s+(\S+)\s*$/){print "$_\n" if $6>=$seqthresh{$1}&&$output{$1}}
 elsif(/^Randomizations:/){last}
 else{print"$_\n"}
 }
#for (keys %matlens) {print "$_\t$seqthresh{$_}\n";}

sub get_seqlens {
    my $filename = shift;
    my @seqlens;
    open FILE, $filename or die $!;
    my $seq;
    while (<FILE>) {
	if ( /^\s*>/ ) {
	    push @seqlens, $seq =~ tr/a-zA-Z// if defined $seq;
	    $seq = '';
	} else {
	    $seq .= $_;
	}
    }
    push @seqlens, $seq =~ tr/a-zA-Z// if defined $seq;  # !
    close FILE;
    return @seqlens;
}

# Perform Fisher's exact test
sub fisher {
    my ($a, $b, $c, $d) = @_;
    my $n = $a + $b + $c + $d;
    my $const_bit = lf($a+$b) + lf($c+$d) + lf($a+$c) + lf($b+$d) - lf($n);
    my $p = 0;
    # would be faster but less accurate to calculate the other tail:
    while ($b >= 0 and $c >= 0) {
        $p += exp( $const_bit - lf($a) - lf($b) - lf($c) - lf($d) );
	--$b; --$c; ++$a; ++$d;
    }
    return $p;
}

# log factorial
sub lf {
    $lookup[$_] = $lookup[$_-1] + log($_) for (@lookup..$_[0]);
    return $lookup[$_[0]];
}

#!/usr/bin/perl
use strict;
use warnings;
my $gtf = shift;
open(my $in, $gtf) or die "Can't open $gtf: $!\n";
my %trans;
my $hasCDS;
my $tID;
while(<$in>){
#	chr	start	end	strand	name	gene.id	isoform.id	cdsStart	cdsEnd	exonCount	exonEnds
#1:	chr1	34554	36081	-	FAM138A	FAM138A	ENST00000417324.1	34553	34553	3	34553,35276,35720,	35174,35481,36081,
  next if $_ =~ /^\#/o;
  chomp $_;
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/,$_);
  next unless($type eq 'transcript' or $type eq 'exon' or $type eq 'CDS');
  my @add_attributes = split(";", $attributes);
  # store ids and additional information in second hash
  my %fields = ();
  foreach my $attr ( @add_attributes ) {
     next unless $attr =~ /^\s*(.+)\s(.+)$/;
     my $c_type  = $1;
     my $c_value = $2;
     $c_value =~ s/\"//go;
     if($c_type  && $c_value){
       $fields{$c_type} = $c_value;
     }
  }
  if($type eq 'transcript'){
    $tID = $fields{transcript_id};
    $trans{$tID}{name} = $fields{gene_name};
    $trans{$tID}{chr} = $chr;
    $trans{$tID}{start} = $start;
    $trans{$tID}{end} = $end;
    $trans{$tID}{strand} = $strand;
    $trans{$tID}{geneid} = $fields{gene_id};
    $trans{$tID}{isoformid} = $fields{transcript_id};
  } elsif($type eq 'exon') {
    $tID = $fields{transcript_id};
    if($strand eq '+'){
      $trans{$tID}{exonStarts} .= $start-1 . ",";
      $trans{$tID}{exonEnds} .= "$end,";
    } else {
      if(defined $trans{$tID}{exonStarts}){
        $trans{$tID}{exonStarts} = $start-1 . "," . $trans{$tID}{exonStarts};
        $trans{$tID}{exonEnds} = "$end," . $trans{$tID}{exonEnds};
      }else{
        $trans{$tID}{exonStarts} = $start-1 . ",";
        $trans{$tID}{exonEnds} = "$end,";
      }
    }
  }elsif($type eq 'CDS') {
    $tID = $fields{transcript_id};
    my $cs = $start - 1;
    if($strand eq '+'){
      $trans{$tID}{cdsStarts} .= $cs-1 . ",";
      $trans{$tID}{cdsEnds} .= "$end,";
    } else {
      if(defined $trans{$tID}{cdsStarts}){
        $trans{$tID}{cdsStarts} = $cs-1 . "," . $trans{$tID}{cdsStarts};
        $trans{$tID}{cdsEnds} = "$end," . $trans{$tID}{cdsEnds};
      }else{
        $trans{$tID}{cdsStarts} = $cs-1 . ",";
        $trans{$tID}{cdsEnds} = "$end,";
      }
    }
  }
}
close $in;

print "chr\tstart\tend\tstrand\tname\tgene.id\tisoform.id\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\n";
foreach my $tID (keys %trans){
  my $tmp = $trans{$tID}{start} - 1;
  my $cs = defined $trans{$tID}{cdsStarts} ? $trans{$tID}{cdsStarts} : $tmp;
  my $ce = defined $trans{$tID}{cdsEnds} ? $trans{$tID}{cdsEnds} : $tmp;
  my $ec = scalar split(/\,/,$trans{$tID}{exonStarts});
  my $cc = scalar split(/\,/,$cs);
  if($cc > 1){
    $cs = $tmp;
    $ce = $tmp;
  }
  $cs =~ s/\,//go;
  $ce =~ s/\,//go;
  $cs = $tmp if($cs < $tmp);
  print "$trans{$tID}{chr}\t$trans{$tID}{start}\t$trans{$tID}{end}\t$trans{$tID}{strand}\t$trans{$tID}{name}\t$trans{$tID}{geneid}\t$trans{$tID}{isoformid}\t$cs\t$ce\t$ec\t$trans{$tID}{exonStarts}\t$trans{$tID}{exonEnds}\n";
}
exit;

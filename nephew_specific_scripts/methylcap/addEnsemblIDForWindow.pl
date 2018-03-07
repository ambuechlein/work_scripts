#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %options;
my $result = GetOptions(\%options,
                        "gene|g=s",
                        "promoter|p=s",
                        "table|t=s",
                       );

warn "Creating Gene Map\n";
my $geneMap = parseMapping($options{gene});
warn "Creating Promoter Map\n";
my $promMap = parseMapping($options{promoter});

warn "Parsing Table\n";
my $of = $options{table};
$of =~ s/(\.txt|\.tsv)$//go;
$of .= ".geneID.tsv";
open(my $tin, $options{table}) or die "Can't open $options{table}: $!\n";
open(my $out, ">$of") or die "Can't open output $of: $!\n";
my $h = <$tin>;
chomp $h;
my @head = split(/\t/,$h);
my $hid = shift @head;
my $print_h = join("\t", "Window","ENSEMBL ID Gene","ENSEMBL ID Promoter",@head);
print $out "$print_h\n";
while(<$tin>){
  chomp $_;
  my @line = split(/\t/,$_);
  my $id = shift @line;
  $id =~ s/_/-/go;
  if(defined ${$geneMap}{$id}){
    foreach my $g (split(/\|/, ${$geneMap}{$id})){
      my $print = join("\t", $id,$g,"-",@line);
      print $out $print, "\n";
    }
  }
  if(defined ${$promMap}{$id}){
    foreach my $g (split(/\|/, ${$promMap}{$id})){
      my $print = join("\t", $id,"-",$g,@line);
      print $out $print, "\n";
    }
  }

  if( (not defined ${$geneMap}{$id}) and (not defined ${$promMap}{$id}) ){
    my $print = join("\t", $id,"-","-",@line);
    print $out $print, "\n";
  }

}
close $tin;
warn "Finished\n";
exit;

sub parseMapping{
  my $file = shift @_;
  open(my $in, $file) or die "Can't open $file: $!\n";
  my %map;
  while(<$in>){
    chomp $_;
    my @line = split(/\t/,$_);
# chr1:9500-10000
    $line[0] =~ /^(chr\d+)\:(\d+)\-(\d+)$/o;
    my $chr = $1;
    my $s = $2;
    my $e = $3;
    $line[0] = "$chr:" . ($s + 1) . "-$e";
    if(not defined $map{$line[0]}){
      $map{$line[0]} = $line[1];
    } else {
      $map{$line[0]} .= "|$line[1]";
    }
  }
  close $in;
  return(\%map);
}

#!/usr/bin/env perl
use strict;
use warnings;
use Excel::Writer::XLSX;
my $out = shift;

# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new( "$out" );
foreach my $file (@ARGV){
  # Add a worksheet
  my $ws = $file;
  $ws =~s/\.DESeq2_Results.annot.tsv//go;
  $ws =~s/\.DESeq2_Results.tsv//go;
  $ws =~ s/\.(tsv|txt)$//go;
  if(length($ws) > 30){
    $ws = substr($ws, 0,30);
  }
  my $worksheet = $workbook->add_worksheet($ws);
  
  # Write a formatted and unformatted string, row and column notation.
  my $row = 0;
  open(my $in, $file) or die "Can't open $file: $!\n";
  while(<$in>){
    chomp $_;
    my @array = split(/\t/,$_);
    my $col = 0;
    foreach my $el (@array){
      $worksheet->write($row, $col, $el);
      $col++;
    }
    $row++;
  }
  close $in;

}
$workbook->close();

exit;

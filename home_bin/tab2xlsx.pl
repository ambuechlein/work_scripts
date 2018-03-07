#!/usr/bin/env perl
use strict;
use warnings;
use Excel::Writer::XLSX;
my $file = shift;
my $out = $file;
$out =~ s/\.(tsv|txt)$//go;
# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new( "$out.xlsx" );

# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Add and define a format
# $format = $workbook->add_format();
# $format->set_bold();
# $format->set_color( 'red' );
# $format->set_align( 'center' );

# Write a formatted and unformatted string, row and column notation.
my $row = 0;
open(my $in, $file) or die "Can't open $file: $!\n";
# $worksheet->write( $. -1, 0, [split] ) while <$in>;
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
$workbook->close();

exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Spreadsheet::WriteExcel;
# Open the tab-delimited file

open (my $tab, $ARGV[0]) or die "$ARGV[0]: $!";
# Create a new Excel workbook
my $workbook  = Spreadsheet::WriteExcel->new($ARGV[1]);
my $worksheet = $workbook->addworksheet();
# Row and column are zero indexed
my $row = 0;
while (<$tab>) {
   chomp $_;
   # Split on single tab
   my @Fld = split('\t', $_);
   my $col = 0;
   foreach my $token (@Fld) {
       $worksheet->write($row, $col, $token);
       $col++;
   }
   $row++;
}
close $tab;
exit;

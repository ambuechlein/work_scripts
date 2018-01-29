#!/usr/bin/env perl
use strict;
use warnings;
my $dir = shift;
opendir(my $in, $dir) or die "Can't open directory: $!\n";
while (my $file = readdir($in)){
# Transcriptome Analysis Run 11 BLAST_vs_NR_Gene_Coverage.png
  if($file =~ /^(\S+)\.geneID.annot(\S+)$/o){
    my $newname = $1.".splitID".$2;
    print "$file | $newname\n";
    rename($file, $newname); 
  }
}
close $in;

exit;

#!/usr/bin/env perl
use strict; 
use warnings;
use Getopt::Long;
my $dir;
GetOptions(
          'directory=s'  => \$dir,
          );

opendir(my $din, $dir) or die "can't open input directory $dir: $!\n";
my $count = 0;
while(my $input = readdir($din)){
  next unless( -f "$dir/$input" and $input =~ /fastq(.(gz|gzip|bz|gzip))?$/o);
# warn "$dir/$input\n"; next;
  open(INPUT, $input =~ /.gz(ip)?$/ ? "zcat $input |" : $input =~ /.bz(ip)?2$/ ? "bzcat $input |" : $input) || die("Open error: $input");
  # Gb = giga base pairs = 1,000,000,000 bp.
  while (<INPUT>) {
      my $line0= $_;
      my $line1=<INPUT>;
      my $line2=<INPUT>;
      my $line3=<INPUT>;
      $count+=length($line1);
  }
  close INPUT;
}
print $count/1000000000, " Gbp\n";
exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %options;
my $result = GetOptions(\%options,
                        "dir|d=s",
                        "tophat|t=s",
                        "sampleorder|s=s",
                        "out|o=s",
                        "jbrowse|j=s",
                       );
my $dir = $options{dir};
my $tophatdir = $options{tophat};
my $map = $options{sampleorder};
my $jb = $options{jbrowse};

open(my $out, ">$options{out}") or die "Can't open output file: $!\n";
print $out <<EOT
#!/bin/csh
setenv PERL5LIB /home/drusch/CPAN/lib/lib:/home/drusch/CPAN/lib/lib/arch:/home/drusch/lib/perl5/lib/perl/5.14.2/:/var/www/jbrowse/src/perl5:/var/www/jbrowse/extlib/lib/perl5:/home/drusch/jcviBin/bioinformatics/utilities/:/home/drusch/jcviBin/frameshiftCorrection/trunk/FrameshiftCorrection/:/nfs/bio/web/tools/jbrowse/jbrowse-1.10.12/src/perl5
rm -r $jb

echo "refseq"
/nfs/bio/web/tools/jbrowse/htdocs/bin/prepare-refseqs.pl --fasta /nfs/labs/nephew/human_databases/hg19/hg19.fasta -out $jb

echo "genes"
/nfs/bio/web/tools/jbrowse/htdocs/bin/flatfile-to-json.pl --gff /nfs/labs/nephew/human_databases/gencode_v19/jbrowse.gencode.v19.annotation.gff --trackLabel "Genes" -out $jb -autocomplete label -type gene
# /nfs/bio/web/tools/jbrowse/htdocs/bin/flatfile-to-json.pl --gff /nfs/labs/nephew/human_databases/gencode_v19/jbrowse.gencode.v19.annotation.gff --trackLabel "Exons" -out $jb -autocomplete label -type exon
EOT
;

open(my $min, $map) or die "can't open $map: $!\n";
my %labels;
my @order;
my %def;
while(<$min>){
  chomp $_;
  my @line = split(/\t/, $_);
  $line[1] =~ s/(\(|\))/_/go; $line[1] =~ s/_+/_/go;
  $labels{$line[0]} = $line[1];
  push(@order, $line[1]);
}
close $min;

print $out "\n", 'echo "Bedgraphs"', "\n";
foreach my $pre (reverse @order){
  my $filen = "$pre.negative.bedgraph";
  my $filep = "$pre.positive.bedgraph";
  warn "$filen does not exist\n" unless(-e "$dir/$filen");
  warn "$filep does not exist\n" unless(-e "$dir/$filep");
  my $color1 = "red";
  my $color2 = "black";
  if($pre =~ /C2D8/go){
    $color1 = "yellow";
    $color2 = "blue";
  }
  print $out "perl /home/drusch/bin/bedGraph2Jbrowse2.pl $jb $filep /nfs/labs/nephew/human_databases/gencode19_with_contaminates/chrSizes.tsv \"$pre positive\" \"$pre positive\" linear local 0 \"$color1\" 30 50\n";
  print $out "perl /home/drusch/bin/bedGraph2Jbrowse2.pl $jb $filen /nfs/labs/nephew/human_databases/gencode19_with_contaminates/chrSizes.tsv \"$pre negative\" \"$pre negative\" linear local 0 \"$color2\" 30 50\n";
}

print $out "\n", 'echo "Junctions"', "\n";
opendir(my $tin, $tophatdir) or die "$!\n";
while(my $file = readdir($tin)){
  next unless($file =~ /.nonredundant\.fasta_thout$/o);
  $file =~ /(\S+)\.fastq\.nonredundant\.fasta_thout/o;
  my $label = "$labels{$1} Junctions";
  my $bg = $tophatdir . '/'. $file . '/junctions.bed';;
  my $color = 'black';

  print $out "/nfs/bio/web/tools/jbrowse/htdocs/bin/flatfile-to-json.pl --bed $bg -out  $jb --trackLabel \"$label\"\n";
}
close $tin;

print $out <<EOT

echo "names"
/nfs/bio/web/tools/jbrowse/htdocs/bin/generate-names.pl -out $jb
EOT
;
close $out;
exit;

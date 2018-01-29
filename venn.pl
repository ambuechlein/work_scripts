#!/usr/bin/env perl
use strict;
use warnings;
use Venn::Chart;

my $asti_b = '/home/abuechle/Desktop/pam_asti_comparison/blast_databases/Asticcacaulis_biprosthecium.faa';
my $asti_e = '/home/abuechle/Desktop/pam_asti_comparison/blast_databases/Asticcacaulis_excentricus.faa';
my $caul = '/home/abuechle/Desktop/pam_asti_comparison/blast_databases/CP001340.faa';

my $a_b = '/home/abuechle/Desktop/pam_asti_comparison/blast_results/Asticcacaulis_biprosthecium_vs_Asticcacaulis_excentricus.txt';
my $a_c = '/home/abuechle/Desktop/pam_asti_comparison/blast_results/Asticcacaulis_biprosthecium_vs_CP001340.txt';
my $b_c = '/home/abuechle/Desktop/pam_asti_comparison/blast_results/Asticcacaulis_excentricus_vs_CP001340.txt';

my @asti_biprosthecium = get_fasta_gene($asti_b);
my @asti_excentricus = get_fasta_gene($asti_e);
my @caulobacter_crescentus = get_fasta_gene($caul);

add_genes($a_b, \@asti_biprosthecium, \@asti_excentricus);
add_genes($a_c, \@asti_biprosthecium, \@caulobacter_crescentus);
add_genes($b_c, \@asti_excentricus, \@caulobacter_crescentus);

my $venn_chart = Venn::Chart->new( 800, 600 ) or die("error : $!");
$venn_chart->set_options( -title => 'Best Blast Hits' );
$venn_chart->set_legends( 'Asticcacaulis biprosthecium', 'Asticcacaulis excentricus', 'Caulobacter crescentus' );
my $gd_venn = $venn_chart->plot( \@asti_biprosthecium,, \@asti_excentricus, \@caulobacter_crescentus );
# Create a Venn diagram image in png, gif and jpeg format
open my $fh_venn, '>', 'VennChart.png' or die("Unable to create png file\n");
binmode $fh_venn;
print {$fh_venn} $gd_venn->png;
close $fh_venn or die('Unable to close file');

# Create an histogram image of Venn diagram (png, gif and jpeg format)
my $gd_histogram = $venn_chart->plot_histogram;
open my $fh_histo, '>', 'VennHistogram.png' or die("Unable to create png file\n");
binmode $fh_histo;
print {$fh_histo} $gd_histogram->png;
close $fh_histo or die('Unable to close file');

# Get data list for each intersection or unique region between the 3 lists
#my @ref_lists = $venn_chart->get_list_regions();
#my $list_number = 1;
#foreach my $ref_region ( @ref_lists ) {
#  print "List $list_number : @{ $ref_region }\n";
#  $list_number++;
#}

exit;

sub get_fasta_gene{
  my $file = shift;
  open(my $in, $file) or die "Can't open $file: $!\n";
  my @genes;
  while(<$in>){
    if($_ =~ /^\>/o){
      $_ =~ s/^\>//o;
      my ($gid) = split(/\s/, $_);
      push(@genes, $gid);
    }
  }
  return @genes;
}

sub add_genes{
  my $file = shift;
  my $array1 = shift;
  my $array2 = shift;
  my %best;
  my %max;
  my %hits;
  open(my $in, $file) or die "Can't open $file: $!\n";
  while(<$in>){
    chomp $_;
    my @F=split /\t/, $_;
    my ($query, $score) = @F[0, 11];

    if (! exists($max{$query}) || $score > $max{$query}) {
      $max{$query} = $score;
      $best{$query} = ();
    }
    if ($score == $max{$query}) {
      push(@{$best{$query}}, $F[1]);
    }
  }
  foreach my $gene (keys %best){
    foreach (@{$best{$gene}}){
      $hits{$_} = 1;
    }
  }
  push(@{$array1}, keys %hits);
  push(@{$array2}, keys %best);
}

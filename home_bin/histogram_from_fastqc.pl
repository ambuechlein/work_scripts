#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use GD::Graph;
use GD::Graph::bars;
use GD::Graph::colour;

my %options;
my $results = GetOptions (\%options,
                          "directory|d=s",
                          "output|o=s",
                         );
my $dir = $options{directory};
opendir(my $din, $dir) or die "Can't open input directory: $!\n";
my $switch = 0;
my %histo;
while(my $sub = readdir($din)){
  next unless(-d $sub and $sub =~ /fastqc$/o);
  opendir(my $sin, "$dir/$sub") or die "Can't open sub directory: $!\n";
  while(my $file = readdir($sin)){
    next unless( $file eq "fastqc_data.txt");
    open(my $in, "$dir/$sub/$file") or die "Couldn't open $file for $sub: $!\n";
    while(<$in>){
      next if($_ =~ /^\#/o);
      chomp $_;
      if($_ =~ /^\>\>Sequence Length Distribution/o){
        $switch = 1;
      } elsif ($_ eq '>>END_MODULE'){
        $switch = 0;
      } elsif ($switch == 1){
        my @line = split(/\t/, $_);
        $histo{$line[0]} += $line[1];
      }
    }
    close $in;
  }
  close $sin;
}
close $din;
open(my $out, ">$options{output}.tsv") or die "can't open output file: $!\n";
my @y; my @x;
print $out "Length\tCount\n";
foreach my $k (sort {$a <=> $b} keys %histo){
  print $out "$k\t$histo{$k}\n";
  push(@y, $histo{$k});
  push(@x, $k);
}

my @data = ([@x], [@y]);
my $graph = GD::Graph::bars->new(800, 600);

GD::Graph::colour::add_colour(cappuccino => [97, 131, 45]);
GD::Graph::colour::add_colour(seaweed => [88, 98, 82]);

$graph->set(
    x_label           => 'Read Length (bp)',
    x_labels_vertical => 1,
    y_label           => 'Count',
    y_tick_number     => 5,
    transparent       => 1,
    fgclr             => 'seaweed',
    dclrs             => [ qw(cappuccino) ],
    accentclr         => 'seaweed',
    shadow_depth      => 1,
    x_label_position  => 1/2,
    bargroup_spacing  => 3,
    title             => 'Read Length Distribution',
) or die $graph->error;

$graph->set_title_font('/home/abuechle/Desktop/CALIBRI.TTF', 30);
$graph->set_x_label_font('/home/abuechle/Desktop/CALIBRI.TTF', 18);
$graph->set_y_label_font('/home/abuechle/Desktop/CALIBRI.TTF', 18);
$graph->set_x_axis_font('/home/abuechle/Desktop/arial.ttf', 10);
$graph->set_y_axis_font('/home/abuechle/Desktop/arial.ttf', 10);
$graph->set_text_clr('seaweed');
my $gd = $graph->plot(\@data) or die $graph->error;

open(IMG, ">$options{output}") or die "Unable to open file: $!";
binmode IMG;
print IMG $gd->png;
exit;

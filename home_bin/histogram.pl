#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use POSIX qw(ceil floor);

my %options = ();
my $results = GetOptions (\%options,
                          'fasta|f=s',
                          'output|o=s');

open(my $in, $options{'fasta'}) or die "can't open input\n";
my @querysize;
while (<$in>){
  next unless $_ =~ /^\>/o;
  #>isotig00001  gene=isogroup00001  length=1905  numContigs=11
  #>TR1|c0_g1_i1 len=288 path=[531:0-287] [-1, 531, -2]
  $_ =~ / len\=(\d+) /g;
  push(@querysize, $1);
}
close $in;

# parse general summary data
&histogram(100, 'Histogram of Sizes of Input Sequences', 'Input Size (bps)', @querysize);


sub histogram{
  my ($bin_width, $title, $xaxis, @list) = @_;
  use GD::Graph;
  use GD::Graph::bars;
  use GD::Graph::colour;
  # This calculates the frequencies for all available bins in the data set
  my %histogram;
  for(@list){
    if((ceil(($_ + 1) / $bin_width) -1) < 40) {
      $histogram{ceil(($_ + 1) / $bin_width) -1}++;
    } else {
      $histogram{40}++;
    }
  }

  my $max;
  my $min;

  # Calculate min and max
  while ( my ($key, $value) = each(%histogram) ){
    $max = $key if !defined($min) || $key > $max;
    $min = $key if !defined($min) || $key < $min;
  }
  my @y;
  for (my $i = $min; $i <= $max; $i++){
    my $bin       = sprintf("% 10d", ($i) * $bin_width);
    my $frequency = $histogram{$i} || 0;
    push(@y, $frequency);
    $frequency = "#" x $frequency;
  }

  my @x;
  my $test = $max*$bin_width > 4000 ? 4000 : $max*$bin_width;
  for(my $i = $min*$bin_width; $i <= $max*$bin_width; $i+=$bin_width){
    push(@x, $i);
  }
  pop(@x);
  push(@x, '4000+');
  my @data = ([@x], [@y]);
  my $graph = GD::Graph::bars->new(800, 600);

  GD::Graph::colour::add_colour(cappuccino => [97, 131, 45]);
  GD::Graph::colour::add_colour(seaweed => [88, 98, 82]);

  $graph->set(
      x_label           => $xaxis,
      x_labels_vertical => 1,
      x_max_value       => 4000,
      y_label           => 'Count',
      y_tick_number     => 10,
      transparent       => 1,
      fgclr             => 'seaweed',
      dclrs             => [ qw(cappuccino) ],
      accentclr         => 'seaweed',
      shadow_depth      => 1,
      x_label_position  => 1/2,
      bargroup_spacing  => 3,
      title             => $title
  ) or die $graph->error;

  $graph->set_title_font('/home/abuechle/Desktop/CALIBRI.TTF', 30);
  $graph->set_x_label_font('/home/abuechle/Desktop/CALIBRI.TTF', 18);
  $graph->set_y_label_font('/home/abuechle/Desktop/CALIBRI.TTF', 18);
  $graph->set_x_axis_font('/home/abuechle/Desktop/arial.ttf', 10);
  $graph->set_y_axis_font('/home/abuechle/Desktop/arial.ttf', 10);
  $graph->set_text_clr('seaweed');
  my $gd = $graph->plot(\@data) or die $graph->error;

  my $basename = "$title.png";
  $basename =~ s/ /_/go;
  open(IMG, ">$options{output}") or die "Unable to open file: $!";
  binmode IMG;
  print IMG $gd->png;
}

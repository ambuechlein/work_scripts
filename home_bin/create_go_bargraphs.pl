#!/usr/bin/env perl

=head1 NAME

create_go_bargraphs.pl - produces bar graphs of the slim categories using count file produced by the map2slim perl script when run using the -c option.

=head1 SYNOPSIS

USAGE: transsummary2yaml.pl 
            --input=/path/to/input_file.txt
            --output=/path/to/output_prefix
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input,-i>
    The input count file produced by the map2slim perl script.

B<--output,-o>
    The output director/prefix for the bar charts produced.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script produces bar graphs of the slim categories. The input is the count file produced by 
the map2slim perl script when run using the -c option.  In theory any goslim may be used, however, 
we use generic goslim maintained by GO.

=head1  INPUT

The input count file produced by the map2slim perl script.  This file is produced by 
running the map2slim script using the -c flag.

=head1  OUTPUT

The output path directory and prefix to use for the output bar graphs.

=head1  CONTACT

    Aaron Buechlein
    abuechle@cgb.indiana.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use SVG::TT::Graph::Bar;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
                          'output|o=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

open(my $in, $options{input}) or logger->logdie("Can't open input file: $!\n");

my @data_b_c;
my @data_b_k;
my @data_b_v;

my @data_c_k;
my @data_c_v;

my @data_m_k;
my @data_m_v;

while(<$in>){
  chomp $_;
  my @line = split(/\t/, $_);
  $line[0] =~ /^GO\:\d+ (.+) \(.+\)$/o;
  my $goname = $1;
  if ($line[4] eq 'biological_process'){
    push(@{$data_b_k[0]}, $goname);
    push(@{$data_b_k[1]}, $line[1]);
  } elsif ($line[4] eq 'cellular_component') {
    push(@{$data_c_k[0]}, $goname);
    push(@{$data_c_k[1]}, $line[1]);
  } elsif ($line[4] eq 'molecular_function') {
    push(@{$data_m_k[0]}, $goname);
    push(@{$data_m_k[1]}, $line[1]);
  } else {
    warn "Woops! Not sure what catgory this is: $line[4]\n";
  }
}

&bar_chart('GO Biological Function', 'biological_function', @data_b_k) if(scalar @data_b_k > 0);
&bar_chart('GO Biological Cellular Location', 'cellular_location', @data_c_k) if(scalar @data_c_k > 0);
&bar_chart('GO Biological Molecular Function', 'molecular_function', @data_m_k) if(scalar @data_m_k > 0);

exit;

sub bar_chart {
  my ($title, $xaxis, @data) = @_;
  use GD::Graph;
  use GD::Graph::bars;
  use GD::Graph::colour;

  my $graph = GD::Graph::bars->new(1024, 768);

  GD::Graph::colour::add_colour(isga_green => [97, 131, 45]);
  GD::Graph::colour::add_colour(grey => [88, 98, 82]);

  $graph->set(
        x_label           => $xaxis,
        x_labels_vertical => 1,
        y_label           => 'Count',
        y_tick_number     => 10,
        transparent       => 1,
        fgclr             => 'grey',
        dclrs             => [ qw(isga_green) ],
        accentclr         => 'grey',
        shadow_depth      => 1,
        x_label_position  => 1/2,
        bar_width         => 30,
        show_values       => 1,
        title             => $title
  ) or die $graph->error;

  $graph->set_title_font('/nfs/bio/db/fonts/DroidSans.ttf', 30);
  $graph->set_x_label_font('/nfs/bio/db/fonts/DroidSans.ttf', 18);
  $graph->set_y_label_font('/nfs/bio/db/fonts/DroidSans.ttf', 18);
  $graph->set_x_axis_font('/nfs/bio/db/fonts/DejaVuSans.ttf', 10);
  $graph->set_y_axis_font('/nfs/bio/db/fonts/DejaVuSans.ttf', 10);
  $graph->set_text_clr('grey');

  my $gd = $graph->plot(\@data) or die $graph->error;
  my $basename = "$title.png";
  $basename =~ s/ /_/go;
  open(IMG, ">$options{output}/$basename") or die( "Unable to open file: $!" );
  binmode IMG;
  print IMG $gd->png;
  close IMG;
}

sub check_parameters {
    my $options = shift;

    unless (-e $options{input}) {
        die("input option not passed or does not exist!");
        exit(1);
    }

    unless (defined $options{output}) {
        die("output option not passed");
        exit(1);
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }
}


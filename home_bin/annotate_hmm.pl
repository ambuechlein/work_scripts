#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Fcntl qw( O_RDONLY );
use MLDBM "DB_File";
use Data::Dumper;

my %args;
my $results = GetOptions (\%args,
                          'input_file|i=s',
                          'output|o=s',
                          'hmm_info|h=s',
#                          'tigr_roles_db_dir|t=s',
#                          'tigrfams_dir|d=s',
                          );

my %hmm_info;
## create the tied hash
tie(%hmm_info, 'MLDBM', $args{hmm_info} ) || die "failed to tie: $!";

open(my $in, $args{input_file}) or die "Can't open input file, $!\n";
<$in>;
<$in>;
my %annotation;
while(<$in>){
  chomp($_);
  my @line = split(/\t/, $_);
  # genomic prediction id
  my $ref_id = $line[8];
  # hmm id
  my $comp_id = $line[10];
  my $total_score = $line[1];
  if ( $total_score >= $hmm_info{$comp_id}{'trusted_cutoff'} ) {
    $annotation{$ref_id}{ 'gene_product_name' } = $hmm_info{$comp_id}{'hmm_com_name'};
    $annotation{$ref_id}{ 'gene_symbol' }       = $hmm_info{$comp_id}{'gene_symbol'};
    $annotation{$ref_id}{ 'EC' }                = $hmm_info{$comp_id}{'ec_num'};
    $annotation{$ref_id}{ 'GO' }                = $hmm_info{$comp_id}{'go'};
    $annotation{$ref_id}{ 'isotype' }           = "HMM::".$hmm_info{$comp_id}{'isotype'};
    $annotation{$ref_id}{ 'TIGR_Role' }         = $hmm_info{$comp_id}{'role_id'};
    $annotation{$ref_id}{ 'hmm_acc' }           = $comp_id;
    $annotation{$ref_id}{ 'trusted_cutoff' }    = $hmm_info{$comp_id}{'trusted_cutoff'};
    $annotation{$ref_id}{ 'total_score' }       = $total_score;

  }

}
close $in;
open(my $out, ">$args{output}") or die "Can't open output file, $!\n";
foreach my $est (sort keys %annotation){
  print $out "$est\t$annotation{$est}{hmm_acc}\t$annotation{$est}{total_score}\t$annotation{$est}{trusted_cutoff}\t@{$annotation{$est}{GO}}\t@{$annotation{$est}{TIGR_Role}}\t$annotation{$est}{isotype}\n";
}

## for my $acc ( keys %info ) {
##     print "$acc\n";
##     
##     for my $key ( keys %{$info{$acc}} ) {
##         if(ref($info{$acc}{$key}) eq 'ARRAY'){
##           print "\t$key\t=>\t@{$info{$acc}{$key}}\n";
##         } else {
##           print "\t$key\t=>\t$info{$acc}{$key}\n";
##         }
##     }
##     print "\n\n";
## }

untie(%hmm_info);

exit(0);


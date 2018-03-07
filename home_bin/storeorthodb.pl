#!/usr/bin/env perl -w 
use Storable qw(nstore_fd retrieve_fd);
use Fcntl qw(:DEFAULT :flock);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s');

#my $orthodbmap = '/nfs/bio/db/OrthoDB/OrthoDB4/tabtext/orthodb4_Vertebrates_tabtext.txt';
my $orthodbmap = $options{input};
open(my $map, $orthodbmap) or die "Can't open $orthodbmap: $!\n";

### Remove the test file in case it exists already ... 
#unlink("OrthoDBGeneMappings.dat");
#unlink("OrthoDBOrthoIDMappings.dat");

my (%genemappings, %orthoidmappings);
my $header = <$map>;
while (my $gene = <$map>){
  chomp $gene;
  my @line = split('\t', $gene);
  $genemappings{$line[2]}{$line[1]} = $gene;
#  push(@{$orthoidmappings{$line[1]}}, $line[2]);
  push(@{$orthoidmappings{$line[1]}}, $gene);
  undef @line;
}
close $map;

$options{input} =~ /OrthoDB5_(\S+)_tabtext/o;
my $output = "OrthoDB".$1;
my $genemappingsfile = $output . "GeneMappings.dat";
#sysopen(GDF, "OrthoDBVertebratesGeneMappings.dat", O_RDWR|O_CREAT, 0666) or die "Can't open OrthoDBGeneMappings.dat: $!\n";
sysopen(GDF, $genemappingsfile, O_RDWR|O_CREAT, 0666) or die "Can't open OrthoDBGeneMappings.dat: $!\n";
flock(GDF, LOCK_EX) or die "Can't lock OrthoDB.dat: $!\n";
nstore_fd(\%genemappings, *GDF) or die "Can't store has to OrthoDB.dat: $!\n";
truncate(GDF, tell(GDF));
close(GDF);

my $idmappings = $output . "OrthoIDMappings.dat";
#sysopen(ODF, "OrthoDBVertebratesOrthoIDMappings.dat", O_RDWR|O_CREAT, 0666) or die "Can't open OrthoDBOrthoIDMappings.dat: $!\n";
sysopen(ODF, $idmappings, O_RDWR|O_CREAT, 0666) or die "Can't open OrthoDBOrthoIDMappings.dat: $!\n";
flock(ODF, LOCK_EX) or die "Can't lock OrthoDB.dat: $!\n";
nstore_fd(\%orthoidmappings, *ODF) or die "Can't store has to OrthoDB.dat: $!\n";
truncate(ODF, tell(ODF));
close(ODF);

# test the stored file
#open(RDF, "OrthoDB.dat") or die "Can't open OrthoDB.dat: $!\n";
#my $rmappings = retrieve_fd(*RDF);
#close(RDF);
#use Data::Dumper;
#print Dumper($rmappings);

exit;

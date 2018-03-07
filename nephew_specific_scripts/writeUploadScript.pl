#!/usr/bin/env perl
use strict;
use warnings;

my $dir = shift;

open(my $out, ">upload2.sh") or die "Can't open upload2.sh: $!\n";
open(my $din, $dir) or die "Can't open $dir: $!\n";
print $out '#!/bin/bash',"\n";
while(<$din>){
  chomp $_;
# #!/bin/bash
# curl -1 -v --disable-epsv --ftp-skip-pasv-ip -u abuechle@indiana.edu:IGUEl93z2aKt --ftp-ssl --ftp-create-dirs --upload-file PTR_Ctl_Exp1_Exp2_5000/PTR_Ctl_Exp1_Exp2_5000.choosePower.tsv ftps://ftp.box.com:990/junco_network_holder/PTR_Ctl_Exp1_Exp2_5000/
#  warn "Skipping $f\n" unless($f =~ /\.(pdf|tsv)$/o);
#  next unless($f =~ /\.(pdf|tsv)$/o);
# spleen_Ctl_Exp1_Exp2_noRep.5000/spleen_Ctl_Exp1_Exp2_noRep.5000.moduleEigengeneClustering.pdf
# spleen_Ctl_Exp1_Exp2_noRep.5000/spleen_Ctl_Exp1_Exp2_noRep.5000GeneLists/spleen_Ctl_Exp1_Exp2_noRep.5000.turquoise_GeneList.tsv
  my $f = $_;
  my $dir = $_;
  $dir =~ s/[^\/]+$//;
#  warn "$dir\n";
#  next if($f =~ /\.RData$/o);
#  next unless($f =~ /(Cytoscape|diffK_table.final.tsv)/);
  print $out "curl -1 -v --disable-epsv --ftp-skip-pasv-ip -u abuechle\@indiana.edu:IGUEl93z2aKt --ftp-ssl --ftp-create-dirs --upload-file $f ftps://ftp.box.com:990/GSF1368_Kim/fastq/$dir\n";
}
close $out;
close $din;
system("chmod 755 upload2.sh");
exit;

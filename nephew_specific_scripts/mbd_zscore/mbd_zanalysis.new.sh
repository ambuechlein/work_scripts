# fasta file containing relevant chromosome sequences (see cmds)
fasta=$1
# regional window size (eg: 25000)
winsize=$2
# regional window step size (eg: 5000)
stepsize=$3
# local bin size (eg: 500) 
binsize=$4
# Sample sheet
samples=$5
# Experimental design
experiment=$6

echo "Converting bam files to density files"
perl bam2density.SGE.pl -d `pwd` -o `pwd`

echo "Computing chromosome lengths"
./chrlen.sh ${fasta}

echo "Binning each density file to regional and local windows"
./binstats.sh ${winsize} ${stepsize} ${binsize}

echo "Exporting regional window files to GFF format"
./bin2gff.sh

echo "Matching local bins with corresponding overlapping regional windows"
./local2regional.sh

echo "Computing Zscores for each local bin with respect to the corresponding overlapping regional windows"
./getZscores.sh

echo "Averaging Zscores across replicates for each sample"
./getZaverages.sh ${samples}

echo "Calculating Z difference (Z delta) between any two samples (Control and Treatment) to be compared"
./getZdiff.sh ${experiment}

echo "Calculating Z ratio from Z difference"
./getZratios.sh

echo "Assigning pvalues to Z ratios"
./getpvalues.sh

echo "Adjusting pvalues for multiple testing by computing FDR"
./getFDR.sh

echo "Done! Final results saved to results folder"

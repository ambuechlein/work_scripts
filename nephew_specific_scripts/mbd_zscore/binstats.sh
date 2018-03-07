#!/bin/bash
# Bin each density file into two sets
# 	- Divided into regional windows, typically of width 25KB with a stepsize of 5KB
# 	- Divided into local bins, typically of width 500 nt


# regional window size (eg: 25000)
winsize=$1
# regional window step size (eg: 5000)
stepsize=$2
# local bin size (eg: 500) 
binsize=$3

# cd to individual chromosome data folder
cd byChr


############ regional window stats ############
# 1. Chromosome
# 2. Bin_ID
# 3. Bin_Start
# 4. Bin_End
# 5. Average_Bin_Coverage
# 6. Standard_Deviation
# 7. Index

## chr1	12500	1	25000	0	0	chr1_1
## chr1	17500	5001	30000	0	0	chr1_2
## chr1	22500	10001	35000	0	0	chr1_3
## chr1	27500	15001	40000	0	0	chr1_4
## chr1	32500	20001	45000	0	0	chr1_5
## chr1	37500	25001	50000	0	0	chr1_6
## chr1	42500	30001	55000	0	0	chr1_7
## chr1	47500	35001	60000	0	0	chr1_8
## chr1	52500	40001	65000	0	0	chr1_9
## chr1	57500	45001	70000	0	0	chr1_10

mkdir -p regional
for c in *.length
do 
	chr=${c%.length}
	for f in ${chr}_*.density
	do 
		fname=${f%.density}
		/home/mnrusimh/Scripts/General/slidingWindow_completeStats.pl --sizes ${c} --input ${f} --chr_col 1 --pos_col 2 --score_col 3 --winsize ${winsize} --stepsize ${stepsize} --output regional/${fname}_binned_regional.txt
	done
done


############ local bin stats ############
# 1. Chromosome
# 2. Bin_ID
# 3. Bin_Start
# 4. Bin_End
# 5. Average_Bin_Coverage
# 6. Index

## chr1	250	1	500	0	chr1_1
## chr1	750	501	1000	0	chr1_2
## chr1	1250	1001	1500	0	chr1_3
## chr1	1750	1501	2000	0	chr1_4
## chr1	2250	2001	2500	0	chr1_5
## chr1	2750	2501	3000	0	chr1_6
## chr1	3250	3001	3500	0	chr1_7
## chr1	3750	3501	4000	0	chr1_8
## chr1	4250	4001	4500	0	chr1_9
## chr1	4750	4501	5000	0	chr1_10

mkdir -p local
for c in *.length
do 
	chr=${c%.length}
	for f in ${chr}_*.density
	do 
		fname=${f%.density}
		/home/mnrusimh/Scripts/General/slidingWindow_stats.pl --sizes ${c} --input ${f} --chr_col 1 --pos_col 2 --score_col 3 --winsize ${binsize} --stepsize ${binsize} --output local/${fname}_binned_local.txt
	done
done
cd ..

#!/bin/bash
<<COMMENT
-----------------------
positional argument:
1st: path/to/SCS file
2nd: path/to/FISH data
3rd: path/to/save/simulate/result
-----------------------
optional arguments:
date: name of directory, use date to distinguish each case
tumorName: name of tumor, default is GBM07
cellNums: total k cell to infer, default is 6
tumorNums: total n tumor samples to deconvolve, default is 3
simuNums: total N data to simulate, default is 20
cellNoise: noise level to add noise to the reference cells
ploidy: 2 or random 
        set the ploidy for clone, default is 2, 
        if the input is random, the simulate the 
        ploidy from uniform distribution in [1,...,8]

COMMENT

python ../MIMOsolver/simulation.py ../data/simulated_$2_integer_CNV.csv\
                                ../data/probe_position.txt\
                                ../simulation\
                                --date $1\
                                --tumorName $2\
                                --cellNums $3\
                                --tumorNums $4\
                                --simuNums $5\
                                --cellNoise $6\
                                --ploidy $7\
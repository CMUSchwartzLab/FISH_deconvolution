#!/bin/bash
<<COMMENT
-----------------------
positional argument:
1st: path/to/interval file
2nd: path/to/data
3rd: name directory to save result, use date to distinguish
-----------------------
optional arguments:
alphaf: regularization for frequency penalty
alpha1: regularization for tree penalty
alpha2: regularization for FISH penalty
offset: index of the data to start from
limit: number of data to run, starting from data[offset]
mask: set to 1
    percentage of genomic interval to be chosen, we use all avaiable genomic loci for now

COMMENT
python ../MIMOsolver/main.py ../data/interval.pkl\
                            ../simulation/$1/GBM07/3\
                            $1\
                            $2\
                            --alphaf $3\
                            --alpha1 $4\
                            --alpha2 $5\
                            --offset $6\
                            --limit  $7\
                            --solver $8\
                            --mask $9\

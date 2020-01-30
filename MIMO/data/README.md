# Data simulation directory

For each of the two tumors GBM07 and GBM33, this directory contains
- three input files (in csv format) for Simulate_SCS.py
- one sample output file (in csv format) from Simulate_SCS.py

The input files are based on the read Single-Cell Sequencing (SCS)
data. GBM07RegionDivide.csv.
Their meanings are as follows:

- GBM33RegionDivide.csv indicate how many cells are in each of three
  regions of each tumor.

- GBM07CopyNumberRate.csv and GBM33CopyNumberRate.csv show the
  probabilities of each copy number in {0,1,3 4,5,6,7,8,9,10} in each
  of three regions of the tumor.

- GBM07PositionRate.csv and GBM33PositionRate.csv show the
  probabilities that each position has a copy number other than 2.

The output files contain the **Simulated** SCS data based on the summary probability distributions of the real SCS
data, so they have the prefix *simulated_* at the beginning of the
file name.

Much more information can be found in
[../LLSolver/README.md](../LLSolver/README.md)

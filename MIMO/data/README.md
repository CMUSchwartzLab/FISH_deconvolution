# Data simulation directory


There there files (in csv format) are based on the read Single-Cell Sequencing (SCS)
data. GBM07RegionDivide.csv.
Their meanings are as follows:

- GBM33RegionDivide.csv indicate how many cells are in each of three
  regions of each tumor.

- GBM07CopyNumberRate.csv and GBM33CopyNumberRate.csv show the
  probabilities of each copy number in {0,1,3 4,5,6,7,8,9,10} in each
  of three regions of the tumor.

- GBM07PositionRate.csv and GBM33PositionRate.csv show the
  probabilities that each position has a copy number other than 2.

We provided the final **Simulated** SCS data based on the summary probability distributions of the real SCS
data, so they have the prefix *simulated_* at the beginning of the file name.

Much more information can be found in our previous work [here](https://github.com/leovam/SCS_deconvolution/tree/master/schwartzlab/data) and [here](https://github.com/leovam/SCS_deconvolution/tree/master/schwartzlab/LLSolver)


We also provide the infomration about the interval and FISH:
- interval.pkl: a file contain three dataframes that indicate the chromosome, begin position and end position of each genomic location in SCS data, as well as such information for FISH probes
- probe_position.txt: a raw file store the position information for FISH data
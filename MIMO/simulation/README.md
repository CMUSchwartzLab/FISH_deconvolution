# simulation directory

This directory and its subdirectories are used to store the simulated Bulk tumor data either from Simulated SCS data or real observed SCS data.

We included 10 replicates for each tumor sample cases in the corresponding subdirectory, please go to each subdirectory for detail.

The one subdirectory named as `1_30_1` means the Bulk data were simulated on the date of `1_30_2020`. This is to distinguish different simulation dataset and we prefer to use the date as the mark. The way of naming it depends on the user. The last digit `1` indicate the case of experiments, in the paper, we test four cases: without noise without ploiy change; with 10% noise without ploidy change; without noise with ploidy change; with 10% noise with ploidy change, we use 1,2,3,4 to distinguish such different case


In each specific date folder, we included small amount of replicates of simulated Bulk data. Each tumor directory contain 3s samples case, each sample case contain 10 replicates for that sample case.

The directory would be like:
```
simulation/DateFolder/simulated_GBM07/3/simulate_data_0
                                       /simulate_data_1
                                     
                                    ...
                     /simulated_GBM33/3/simulate_data_0
                                       /simulate_data_1
                                     
                                    ...                               
```



While we cannot redistribute human subject data, we still used created `GBM07` and `GBM33` with user's reference. Please keep in mind the simulated data here is not from real SCS data but from simulated SCS data.

Much more information can be found in
[../MIMOsolver/README.md](../MIMOsolver/README.md).
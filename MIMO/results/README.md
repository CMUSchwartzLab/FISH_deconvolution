# results directory

This directory and its subdirectories are used to store the results of copy mumber Deconvolution from the bulk data with information from FISH

The one subdirectory named as `1_30_1_1.0` means the deconvolution model are test on `1_30_2020`. If the user simulated the data and test on it at in the same day, this date would also infer which simulated data are tested on, however, the user should keep record on what simulated data are tested on and use the date here only for the different test.

More information on the date directory and output files:
- **1_30_1_1.0**: 1_30 indicate the date, 1 indicate the case of experiments, in the paper, we test four cases: without noise without ploiy change; with 10% noise without ploidy change; without noise with ploidy change; with 10% noise with ploidy change, we use 1,2,3,4 to distinguish such different case. 1.0 means that we use all the genomic loci as reference
- **sample_output_0_0.2_0.2_0.2.txt**: a raw output of the inference, 0 indicate the index of the output, three 0.2 indicate the value for alpha_f, alpha_1 and alpha_2 respectively. In the txt file are different sections that are seperated by ```---```, from the top to the bottom are the raw output of copy number (9934 *6, 9934 indicates the size of genomic loci, 6 indicates the number of cell components); fractions (6 * 3, 6 indicates the number of cell components and 3 indicates the number of tumor samples); ploidy (an array of 6, each element is the ploidy of each cell component); phylogenetic tree (if there) (13 * 13, first 6 columns are inferred cell components, second 6 the observed cell components and last column the root)
- **stastics_output_0_0.2_0.2_0.2.txt**: statistic calculation for performance for each case

Much more information can be found in
[../MIMOsolver/README.md](../MIMOsolver/README.md) and [../simulation/README.md](../simulation/README.md)
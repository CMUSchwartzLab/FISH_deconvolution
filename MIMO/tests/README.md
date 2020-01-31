# Test Script Directory

This directory contains the scripts that call different .py in **MIMOsolver**. We also provide corresponding calling script in order to save the input argumemnts. When running the tests for the manuscript, such scripts were submitted to a compute farm that accepts bash scripts. 

For example
simulate_command.sh provide a structure of argumemnts to call simulation.py, following the instruction, a user can call simulation.py by running
```bash
bash simulate_command.sh 1_30_1 GBM07 6 3 10 0 2
```
in terminal, or can specify the argument in run_simulate_command.sh and run
```bash 
bash run_simulate_command.sh
```

Much more information can be found in
[../MIMOsolver/README.md](../MIMOsolver/README.md).

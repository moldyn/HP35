# Robust Density-Based Clustering

In this directory you find the two final microstates trajectories used in Nagel et al. 2023 and the script `perform_clustering` to reproduce them.

## Workflow to Create Microstates
Here is a brief summarize of how the two files were generated.
- `hp35.dihs.res3-33.shifted.gaussian10f_microstates_pcs4_p153`
    1. Estimate free energy with default parameters
    1. Apply free energy screening
    1. Apply minimal population cut-off with $P_\text{min}=153$
    1. Create microstate trajectory
- `hp35.mindist2.gaussian10f_microstates_pcs5_p153`
    1. Estimate free energy with default parameters
    1. Apply free energy screening
    1. Apply minimal population cut-off with $P_\text{min}=153$
    1. Create microstate trajectory

For a more detailed description of the overall workflow, please take a look at the publication mentioned above.
For more details on each step, please take a look at the [documentation](https://moldyn.github.io/Clustering/).

## Reproduce the Results
To reproduce the resulting trajectory you simply have to run
```bash
# this will download and compile all needed files
bash perform_clustering -c 1
```
This creates the directory `create_clustering_nagel23` with the two projections including all intermediate steps. If you need more information on the executed commands, you can run the script in the verbose mode `-v`.

# Most Probable Path Lumping (MPP)

In this directory you find the two final macrostates trajectories used in Nagel et al. 2023 and the script `perform_mpp` to reproduce these results.

## Workflow to Create Macrostates
Here we want to briefly summarize how the files were generated. For a more in depth description, please take a look at the publication.

1. Applying most probable path lumping
    - In contrast to PCCA, MPP depends on the free energy. Therefore, we included the execution to the [perform_clustering](../CLUSTERING/perform_clustering) script.
    - In general one obtains the needed files via `clustering mpp -s microstates -l lagtime -D fe -v` where `clustering` can be installed via `conda install moldyn-clustering -c conda-forge`, `microstates` is the microstate file, `fe` the file containing the free energy and the `lagtime` is given in frames.
    - This results in the linkage matrix file, called `*_transitions.dat` 
1. Automated branch detection
    - In Nagel et al. 2023 we suggested an improved version of MPP, which allows to automatically identify branches and perform an optimal lumping.
    - ```python process_mpp.py --linkage linkagemat --tlag lagtime --state-traj microstates --cut-params minpop minq --fraction-of-native-contacts q_of_t --hide-labels```

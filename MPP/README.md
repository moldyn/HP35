# Most Probable Path Lumping (MPP)

In this directory you find the two final macrostates trajectories used in Nagel et al. 2023 and the script `perform_mpp` to reproduce them.

## Workflow to Create Macrostates
Included macrostate trajectories:
1. Dihedral-based macrostates: `hp35.dihs.res3-33.shifted.gaussian10f_microstates_pcs4_p153.mpp50_transitions.dat.renamed_by_q.pop0.001_qmin0.50.macrotraj`
    - Using a lagtime $\tau=10\text{ns}$
    - Using a required minimal population of $0.1%$ and a required minimal metastability of $0.5$.
2. Contact-based macrostates: `hp35.mindists2.gaussian10f_microstates_pcs5_p153.mpp50_transitions.dat.renamed_by_q.pop0.005_qmin0.50.macrotraj`
    - Using a lagtime $\tau=10\text{ns}$
    - Using a required minimal population of $0.5%$ and a required minimal metastability of $0.5$.
3. Contact-based macrostates: `hp35.mindists2.gaussian10f_microstates_pcs5_p153.mpp50_transitions.dat.renamed_by_q.pop0.005_qmin0.50.macrotraj_lumped13`
    - Same as above, but states state 13 was merged into state 12, so $12, 13 \to 12$.

Here is a short summary of how the macrostate trajectories were generated. For a more detailed description of the overall workflow, please take a look at the publication above.

1. Applying most probable path lumping
    - In contrast to PCCA, MPP depends on the free energy. Therefore, we included the execution in the [perform_clustering](../CLUSTERING/perform_clustering) script.
    - In general one obtains the needed files via
      ```
      clustering mpp \  # install via `conda install moldyn-clustering -c conda-forge`
          -s microstates \  # path to microstate trajectory 
          -l lagtime \  # lagtime to estimate Markov model [frames]
          -D fe \  # time series of free energy estimation
          -v  # verbose mode
      ```
    - This results in the linkage matrix file, called `*_transitions.dat` 
2. Automated branch detection
    - In Nagel et al. 2023 we suggest an improved version of MPP, which allows to automatically identify branches and perform an optimal lumping.
      ```bash
      python process_mpp.py \
          --linkage linkagemat \  # path to linkage matrix obtained via mpp
          --tlag lagtime \  # same lagtime as used for mpp [frames]
          --state-traj microstates \  # path to microstate trajectory
          --cut-params minpop minq \  # minimum population and metastability both [0, 1]
          --fraction-of-native-contacts q_of_t \  # time series of fraction of native contacts
          --hide-labels  # hide x-ticklabels
      ```
    - The time evolution of the fraction of native contacts q is available in [create_macrostate_nagel23/hp35.mindists2.gaussian10f.q](create_macrostate_nagel23/hp35.mindists2.gaussian10f.q)
    - Please ensure that all dependencies are installed by creating a new conda environment or using venv.
      ```bash
      conda create -n hp35 -c conda-forge python msmhelper tqdm prettypyplot scipy click numpy
      conda activate hp35
      ```

## Reproduce the Results
To reproduce these results you simply have to run
```bash
# this will create automatically a new Python environment
bash perform_mpp -c 1
```
This creates the dendrograms and the macrostate trajectories in the directory `create_macrostate_nagel23` . If you need more information on the executed commands, please run the script in the verbose mode `-v`.

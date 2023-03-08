# Principal Components

In this directory you find the two final principal component projections used in Nagel et al. 2023 and the script `perform_pca` to reproduce them.

## Workflow to Create PCs
Here is a brief summary of how the two files were generated. For a more detailed description, please take a look at the publication above.

- `hp35.dihs.res3-33.shifted.gaussian10f.proj.1-4`
    1. Discarding angles of residue 2 and 34, $\phi_2, \psi_2, \phi_{34}, \psi_{34}$
    1. Apply Gaussian filtering with $\sigma=2\text{ns}$
    1. Perform PCA, here we do not need using dPCA+ because the dihedrals are already maximal gap shifted
    1. Extract first 4 PCs
- `hp35.mindists2.gaussian10f.proj.1-5`
    1. Apply Gaussian filtering with $\sigma=2\text{ns}$
    1. Perform conPCA
    1. Extract first 5 PCs


## Reproduce the Results
To reproduce these results you simply have to run
```bash
# this will download and compile all needed files
bash perform_pca -c 1
```
This creates the directory `create_pca_nagel23` with the two PCs projections, including all intermediate steps.

# Selecting Features for Markov Modeling: A Case Study on HP35
This repository provides all scripts and intermediate steps to generate the
reproduce the analysis of Nagel et al. 2023. If the provided scripts/files are
used, please cite:
> D. Nagel, S. Sartore, and G. Stock,  
> *Selecting Features for Markov Modeling: A Case Study on HP35*,  
> J. Chem. Theory Comput., submitted;  
> doi: [xx.xxxx/x.xxxxxxx](https://aip.scitation.org/doi/xx.xxxx/x.xxxxxxx)

## Getting Started
To download all included submodules, please clone this repository with
```bash
git clone --recurse-submodules git@github.com:moldyn/HP35.git
cd HP35
```
## Create States
### Features: Backbone Dihedral Angles and Minimal Contact Distances
In the directory `HP35-DESRES` you can find

1. `hp35.dihs`: backbone dihedral angels given [degrees]
1. `hp35.dihs.shifted`: maximum-gap shifted backbone dihedral angels [rad]
1. `hp35.crystaldists`: the atom distances of all contacts occurring in the crystal structure 2f4k [nm]
1. `hp35.mindists`: all minimal distances occurring more frequently than 30% [nm]
1. `hp35.mindists2`: improved distances definition with all atom pairwise distances occurring more frequently than 30% [nm]

For more details take a look at the repository [HP35-DESRES](https://github.com/moldyn/HP35-DESRES). For legal reasons, the data cannot be integrated directly. Please note the attached [license](https://github.com/moldyn/HP35-DESRES/blob/main/LICENSE.txt) when using.

### Principal Components
In the directory [PCA](PCA) you can find the resulting principal component projections

1. `hp35.dihs.res3-33.shifted.gaussian10f.proj.1-4`
1. `hp35.mindists2.gaussian10f.proj.1-5`

Furthermore you can find script to reproduce them. For more information please take a look at the [README](PCA).

### Microstate Trajectories
In the directory [CLUSTERING](CLUSTERING) you can find the resulting microstate trajectories.
1. `hp35.dihs.res3-33.shifted.gaussian10f_microstates_pcs4_p153`
1. `hp35.mindist2.gaussian10f_microstates_pcs5_p153`
Furthermore you can find script to reproduce them. For more information please take a look at the [README](CLUSTERING).

## Markov State Analysis

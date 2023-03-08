# Selecting Features for Markov Modeling: A Case Study on HP35
This repository provides all the scripts and intermediate steps to
reproduce the analysis of Nagel et al. 2023. If the provided scripts/files are
used, please cite:
> D. Nagel, S. Sartore, and G. Stock,  
> *Selecting Features for Markov Modeling: A Case Study on HP35*,  
> J. Chem. Theory Comput., submitted  
> doi: [10.48550/ARXIV.2303.03814](https://doi.org/10.48550/ARXIV.2303.03814)

## Getting Started
To download all the included submodules, please clone this repository with
```bash
git clone --recurse-submodules git@github.com:moldyn/HP35.git
cd HP35
```
## Create States
### Features: Backbone Dihedral Angles and Minimal Contact Distances
In the directory `HP35-DESRES` you can find

1. `hp35.dihs`: backbone dihedral angles [degrees]
1. `hp35.dihs.shifted`: maximum-gap shifted backbone dihedral angles [rad]
1. `hp35.crystaldists`: the atom distances of all contacts occurring in the crystal structure 2f4k [nm]
1. `hp35.mindists`: all minimal distances occurring more frequently than 30% of the time [nm]
1. `hp35.mindists2`: improved contact distances definition with all atom pairwise distances occurring more frequently than 30% of the time [nm]

For more details, take a look at the repository
[HP35-DESRES](https://github.com/moldyn/HP35-DESRES). For legal reasons, the
data cannot be integrated directly. Please note the attached
[license](https://github.com/moldyn/HP35-DESRES/blob/main/LICENSE.txt) when
using it.

### Principal Components
In the directory [PCA](PCA) you can find the resulting principal component projections

1. `hp35.dihs.res3-33.shifted.gaussian10f.proj.1-4`
1. `hp35.mindists2.gaussian10f.proj.1-5`
In the same directory you can also find a script to reproduce them. For more information please
take a look at the [README](PCA).

### Microstate Trajectories
In the directory [CLUSTERING](CLUSTERING) you can find the resulting microstate trajectories.
1. `hp35.dihs.res3-33.shifted.gaussian10f_microstates_pcs4_p153`
1. `hp35.mindist2.gaussian10f_microstates_pcs5_p153`

In the same directory you can also find a script to reproduce them. For more information please
take a look at the [README](CLUSTERING).


### Macrostate Trajectories
In the directory [MPP](MPP) you can find the resulting macrostate trajectories.
1. `hp35.dihs.res3-33.shifted.gaussian10f_microstates_pcs4_p153.mpp50_transitions.dat.renamed_by_q.pop0.001_qmin0.50.macrotraj`
1. `hp35.mindists2.gaussian10f_microstates_pcs5_p153.mpp50_transitions.dat.renamed_by_q.pop0.005_qmin0.50.macrotraj_lumped13`

In the same directory you can also find a script to reproduce them. For more information please
take a look at the [README](MPP).


## Markov State Analysis
In the directory [MSM](MSM) you can find a description of how to reproduce the
Markov state model analysis of Nagel et al. 2023. All the analysis is based on
the Python package [msmhelper](https://moldyn.github.io/msmhelper/).

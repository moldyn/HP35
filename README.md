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
## Intermediate Steps
### Features: Backbone Dihedral Angles and Minimal Contact Distances
In the directory `HP35-DESRES` you can find
1. `hp35.dihs`: backbone dihedral angels given [degrees]
1. `hp35.dihs.shifted`: maximum-gap shifted backbone dihedral angels [rad]
1. `hp35.crystaldists`: the atom distances of all contacts occurring in the crystal structure 2f4k [nm]
1. `hp35.mindists`: all minimal distances occurring more frequently than 30% [nm]
1. `hp35.mindists2`: improved distances definition with all atom pairwise distances occurring more frequently than 30% [nm]
for more details take a look at the [README](HP35-DESRES/README.md).

### Principal Components

### Microstate Trajectories


## Reproducing the Results
### Getting Started
Simply clone this repository with
```bash
git clone --recurse-submodules git@github.com:moldyn/HP35.git
cd HP35
```
### Reproducing the Microstate Trajectories
Reproducing the published microstate trajectories can be achieved with
```bash
cd clustering && bash robust_clustering -c1
```

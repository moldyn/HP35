[![wemake-python-styleguide](https://img.shields.io/badge/style-wemake-000000.svg)](https://github.com/wemake-services/wemake-python-styleguide)
# HP35-Benchmark
Benchmarking dimensionality reduction methods and clustering on [HP35-DESRES](https://github.com/moldyn/HP35-DESRES).

This is package automatizes the comparison of  different clustring and dimensionality reduction methods and its effects on the resulting Markov state models.

If the provided scripts are used, please cite:
- D. Nagel, A. Weber, B. Lickert and G. Stock, *Dynamical coring of Markov state models*, J. Chem. Phys., 150, 094111, 2019; DOI: [10.1063/1.5081767](https://aip.scitation.org/doi/10.1063/1.5081767)

## Getting ready
Simply clone this repository with
```bash
git clone --recurse-submodules git://github.com/moldyn/HP35-Benchmark.git
cd HP35-Benchmark
```

Reproducing the published state trajectory can be achieved with
```bash
cd clustering && bash robust_clustering
```

## Add Own Routine
For an example take a look at [robust_clustering](https://github.com/moldyn/HP35-Benchmark/blob/master/clustering/robust_clustering)

Structural Variation Evolutionary Tree Analysis
=============

### Installation
To build the C++ simulation code, change into the code directory and type make:
```shell
> cd code
> make
```

### Dependencies
* cmake is required
* plot-cns.R requires R/bioconductor packages `copynumber`, `reshape` and `tools`
* plot-trees-all.R requires R/bioconductor packages `ggtree`, `ape` and `tools`

### Running
You may use the provided scripts to run the programs.

```shell
# Simulating mutations on a coalescence tree
> bash run-sveta.sh
# Build a tree from copy number profile with maximum likelihood method
> bash run-svtreeml.sh
# Build a tree from copy number profile with MCMC method
> bash run-svtreemcmc.sh
```

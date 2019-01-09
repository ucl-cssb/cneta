Structural Variation Evolutionary Tree Analysis
=============

### Installation
To build the C++ simulation code, change into the code directory and type make:
```shell
> cd code
> make
```
  
### Dependencies
plot-cns.R requires R/bioconductor packages `copynumber`, `reshape` and `tools`  
plot-trees.R requires R/bioconductor packages `ggtree`, `ape` and `tools`

### Running
The bash scripts show to run the simulations. Essentially all parameters are passed on the command line and must be in the correct order.
dir is output directory, Nr is the number of regions and Ns is the number of tumours to simulate  
r1, r2, r3, r4, r5 are the rates of duplication, deletion, chromosome gain, chromosome loss, wgd
s1, s2 are the mean of the duplication and deletion size distributions

```shell
> code/sveta <dir> <Nr> <Ns> <r1> <r2> <r3> <r4> <r5> <s1> <s2>
> Rscript ana/plot-trees.R <dir>
> Rscript ana/plot-cns.R <dir>
```
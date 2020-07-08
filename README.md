Structural Variation Evolutionary Tree Analysis
=============

# Introduction
This is a set of programs to simulate and build phylogenetic trees from copy number profiles caused by chromosomal alteration events and structural variations (SVs).
Currently, five types of events are considered, including segment duplication, segment deletion, chromosomal gain, chromosomal loss, and whole genome doubling.

The tree building programs take as input the allele-specific or total copy numbers called from mulitiple samples of a patient.

There are mainly 3 programs:
* sveta: simulating SVs along a phylogenetic (coalescence) tree
* svtreeml: building phylogenetic trees from copy numbers with maximum likelihood approach
* svtreemcmc: building phylogenetic trees from copy numbers with Bayesian MCMC approach


# Installation
This package is mostly written in C++. There are a few scripts written in R and Python, for plotting and text processing.

## Dependencies

* Required C/C++ libraries
  * CMake is required for BFGS optimization
  * C libaries: gsl, boost (version >= 1.42)

* Required R libraries
  * plot-cns.R: `copynumber`, `reshape`, `tools`, `tidyr`, `dplyr`, `purrr`
  * plot-trees-all.R: `ggtree`, `ape`, `tools`, `ggplot2`
  * tree_comparison.R: `phangorn`

* Required Python libraries
  * newick2elist.py: networkx

### How to install CMake

Get the latest “Unix/Linux Source” *.tar.gz file.
```
HOME=~/local
mkdir -p $HOME
tar -xf cmake*.tar.gz
cd cmake*
HOME=~/local
./configure --prefix=$HOME
make
make install
```

### How to install required R libraries
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("copynumber")
BiocManager::install("ggtree")

install.packages("reshape")
install.packages("ape")

# Check installed packages
installed.packages()[, c("Package", "LibPath")]

```

## Building C++
OpenMP is used to accelerate tree search in svtreeml.
To turned off OpenMP, please set "omp =" in makefile.

To build the C++ code, change into the code directory and type make:
```shell
> cd code
> make
```


## Running
You may use the provided bash scripts to run the programs.

```shell
# Simulating mutations on a coalescence tree
> bash run-sveta.sh
# Build a tree from copy number profile with maximum likelihood method
> bash run-svtreeml.sh
# Build a tree from copy number profile with MCMC method
> bash run-svtreemcmc.sh
```


# Simulation with sveta
There are three Markov models of evolution for the copy number profiles:
* 0: Mk model (extension of JC69 model)
* 1: model of total copy number
* 2: model of allele-specific copy number
<!-- * 3: model of independent Markov chains (when WGD and chromosome gain/loss are incorporated) -->

There are two ways of simulating mutations along a tree:
1. simulating waiting times along a branch (default)
2. simulating sequences at the end of a branch
Note that when simulating waiting times, only allele-specific model is allowed. If model of total copy number is specified, it will automatically be converted to model of allele-specific copy number.


Please see run-sveta.sh to learn how to set different parameters

## Input
* --epop Ne: Ne is used to scale the tree height so that the branch length is measured in the unit of year, by multipling each branch length with Ne.

* --tiff delta_t: On the initial tree, the tip nodes have the same time. This parameter can be used to introduce different times at the tip nodes. The terminal branches are increased by random multiples of delta_t. The maximum multiple is the number of samples.

* --cn_max cn_max: The maximum copy number allowed in the program depends on the heap space.


## Output
* *-cn.txt.gz: The total copy number for each site on each sample
* *-rel-times.txt: The sampling time of tip nodes
* *-allele-cn.txt.gz: The allele-specific copy number for each site on each sample
* *-info.txt: The time of each node and the total number of mutations simulated on each branch of the tree, grouped by lineages of tip nodes.
* *-mut.txt: The list of simulated mutations on each branch of the tree.
* *-tree.txt: The simulated tree in tab-delimited format
* *-tree.nex: The simulated tree in NEWICK format, with branch length reprenting calender time
* *-tree-nmut.nex: The simulated tree in NEWICK format, with branch length reprenting number of mutations

File *-cn.txt.gz can serve as the input to a tree building program that used integer total copy number.

File *-allele-cn.txt.gz can serve as the input to a tree building program that used integer allele-specific copy number.

File *-rel-times.txt can provide the timing information of tip nodes to allow etimation of divergence time and mutation rates.

Files *-tree.* provide the real tree, which can be used for measuring the accuracy of tree building programs.

File *-info.txt and *-mut.txt can be used to map mutations onto the tree.
The columns in *-mut.txt: sample_ID, edge_ID, muttype_ID, mut_btime, mut_etime, chr_haplotype, chr, seg_ID
For chromosome gain/loss, seg_ID is assigned to -1.
For whole genome doubling, chr is assigned to 0 and seg_ID is assigned to -1.



# Tree building with ML
<!-- ## How to prepare MP trees -->
## Input
The initial trees for tree searching can be obtained by maximum parsimony methods.

There are 4 running modes in svtreeml.
* mode 0: building maximum likelihood tree from input copy numbers
* mode 1: a simple comprehensive test on a simulated tree
* mode 2: computing likelihood given a tree and its parameters (branch length and mutation rates)
* mode 3: computing maximum likelihood tree given a tree topolgy
The last three modes can be used to validate the computation of likelihood.

There are 3 tree searching method:
* exhaustive search (feasible for trees with fewer than 7 samples)
* hill climbing
* genetic algorithm (may be slow)

Please see run-svtreeml.sh to learn how to set different parameters

## Output


## Model of evolution
For tree recontruction with svtreeml on data with chromosome gain/loss and WGD, model 3 (model of independent Markov chains) should be used. This model is different from the model used for simulating the data, which is usually model 2 (model of allele-specific copy number) in sveta.


# Tree building with MCMC

## Input
* (Required) A file containing copy numbers for all the samples, including the normal sample (*-cn.txt.gz or *-allele-cn.txt.gz)
* (Optional) A file containing the timing information of tip nodes (*-rel-times.txt)
* (Optional) A configuration file which sets most input parameters (mcmc.cfg)

Please see run-svtreemcmc.sh to learn how to set different parameters

There are two runing modes depending on whether a reference tree is provided.
With a reference tree, the tree topolgy is fixed.


## Output
There are two output files in a format similar to that of MrBayes:
* *.p, which records the traces of parameters
* *.t, which records the sampled trees

<!-- ## How to analyze the results of svtreemcmc -->
The output can be analyzed by [RWTY](https://github.com/danlwarren/RWTY). Please see script ana/check_convergence.R for reference.

*.p can be imported into [Tracer](https://beast.community/tracer) to check the convergence of the chains.

*.t can be analyzed by [TreeAnnotator](https://beast.community/treeannotator) to get a summary tree (maximum credibility tree).

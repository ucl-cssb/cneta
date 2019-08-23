Structural Variation Evolutionary Tree Analysis
=============

# Introduction
This is a set of program to build phylogenetic trees from mutational profiles caused by structural variations (SVs).
Currently, five types of SV events are considered, including segment duplications, segment deletions, chromosomal gain, chromosomal loss, and whole genome doubling.

The tree building programs take as input the total copy numbers called from patient samples.



There are 3 major programs:
* sveta: simulating SVs along a phylogenetic (coalescence) tree
* svtreeml: building phylogenetic trees from total copy numbers with maximum likelihood approach
* svtreemcmc: building phylogenetic trees from total copy numbers with MCMC approach


# Installation

## Dependencies

This package is mostly written in C++. There are a few scripts written in R and Python, for plotting and text processing.

* Required C/C++ libraries
  * CMake is required for BFGS optimization
  * C libaries: gsl, boost (version >= 1.42)

* Required R libraries
  * plot-cns.R: `copynumber`, `reshape` and `tools`
  * plot-trees-all.R: `ggtree`, `ape` and `tools`

* Required Python libraries
  *

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

### How to required R libraries
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


## Building
To build the C++ code, change into the code directory and type make:
```shell
> cd code
> make
```

## Running
You may use the provided scripts to run the programs.

```shell
# Simulating mutations on a coalescence tree
> bash run-sveta.sh
# Build a tree from copy number profile with maximum likelihood method
> bash run-svtreeml.sh
# Build a tree from copy number profile with MCMC method
> bash run-svtreemcmc.sh
```

### Plotting
bin_locations_4401.Rdata constains the position of 4401 bins in the human reference genome.
The bins are obtained by non-overlapping sliding windows of size 500,000 bp along each chromosme.

# Simulation with sveta

## Input

* --epop Ne: The initial coalescence tree has a expected tree height smaller than 2. Ne can be used to scale the tree height by Ne, by multipling each branch length with Ne.

* --tiff delta_t: On the initial tree, the tip nodes have the same time. This parameter can be used to introduce different times at the tip nodes. The terminal branches are increased by random multiples of delta_t. The maximum multiple is the number of samples. The terminal branch lengths in the initial coalescence tree are usually very small, with expected length being 1/C(Ns,2) 

## Output
* *-cn.txt.gz: The total copy number for each site on each sample
* *-rel-times.txt: The sampling time of tip nodes
* *-allele-cn.txt.gz: The allele-specific copy number for each site on each sample
* *-info.txt: The time of each node and the mutations simulated on each branch of the tree, grouped by lineages of tip nodes.
* *-mut.txt: The list of simulated mutations on each branch of the tree.
* *-tree.txt: The simulated tree in tab-delimited format
* *-tree.nex: The simulated tree in NEWICK format, with branch length reprenting calender time
* *-tree-nmut.nex: The simulated tree in NEWICK format, with branch length reprenting number of mutations

File *-cn.txt.gz can serve as the input to a tree building program that used total copy number.

File *-allele-cn.txt.gz can serve as the input to a tree building program that used allele-specific copy number.

File *-rel-times.txt can provide the timing information of tip nodes to allow etimation of divergence time and mutation rates.

Files *-tree.* provide the real tree, which can be used for measuring the accuracy of tree building programs.

File *-info.txt and *-mut.txt can be used to map mutations onto the tree?



# Note
## Model of evolution
There are 3 models of evolution:
* Mk model
* model of total copy number
* model of allele-specific copy number


## How to prepare MP trees

## Different options in svtreeml
There are 4 running modes in svtreeml.
* mode 0: building maximum likelihood tree from input copy numbers
* mode 1: a simple comprehensive test on a simulated tree
* mode 2: computing likelihood given a tree and its parameters (branch length and mutation rates)
* mode 3: computing maximum likelihood tree given a tree topolgy

The last three modes can be used to validate the computation of likelihood.


There are 3 tree searching method:
* exhaustive search (only feasible for trees with fewer than 7 samples)
* hill climbing
* genetic algorithm (may be slow)





## How to analyze the results of svtreemcmc
1.

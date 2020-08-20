Structural Variation Evolutionary Tree Analysis
=============

# Introduction
This is a set of programs to simulate and build phylogenetic trees from copy number profiles caused by chromosomal alteration events and structural variations (SVs).
Currently, five types of events are considered, including segment duplication, segment deletion, chromosomal gain, chromosomal loss, and whole genome doubling (WGD).

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
SVETA simulates structural variations that alter copy numbers along a coalescence tree of multiple samples.
The coalescence tree can be either normal or exponential growing.
Each node in the tree corresponds to the genome of a sample, which is represented by a consecutive set of pre-specified sites.
Each site is considered as a segment of unknown size.
These sites can be seen as the segments obtained by segmentation methods when calling copy numbers from real data.
The tip dates can be adjusted to reflect different sampling times.

The program generates a set of files, which record the simulated allele-specific and total copy number profiles, tree topology, and mutations along the branches respectively.
The simulated copy number profiles and/or tip timing information served as input for the tree building methods.


The procedure of simulations is as follows:
1. Generate a random coalescence tree (with exponential growth) which represents the genealogy of tumour samples. Available trees can also be given as input, such as trees generated by [CoalEvol](https://code.google.com/archive/p/coalevol/). The tree can be constrained by a specified patient age and the tip nodes can have different times.
		1) Add a branch of length zero which leads to a new top node representing a normal genome.
		2) (optional) Assign random times (in year) to the tips to create different sampling times, or the times can be specified by parameters (TODO).
		3) (optional) Rescale the tree by the specified patient age at first sampling time so that the tree height is no larger than the patient age at last sampling time.
2. Simulate copy number profiles on the tree by using exponential waiting times and the jump chain (copy number substitution model).
		1) Generate the copy number profiles for the root (normal diploid genome).
		2) Simulate mutational events along each branch.
3. Output result files


There are three Markov models of evolution for the copy number profiles:
* model 0: Mk model (extension of JC69 model)
* model 1: bounded model of total copy number
* model 2: bounded model of allele-specific copy number
<!-- * 3: model of independent Markov chains (when WGD and chromosome gain/loss are incorporated) -->

There are two ways of simulating mutations along a tree:
1. Simulating waiting times along a branch (default). A random waiting time is generated from the exponential distribution with mean $1/r$. Here, $r$ is the total mutation rate across the genome, obtained by adding up the duplication and deletion rates across all sites in the genome, chromosomal gain or loss rates across all chromosomes and whole genome doubling rate. When a mutation is generated, its type is randomly chosen based on the relative rates of different types of mutational events.
2. Simulating sequences at the end of a branch. This is appropriate when the mutational events are not of interest. In this way, we can quickly simulate a much larger number of segments with a larger range of mutation rates.

Note that when simulating waiting times, only allele-specific model is allowed. If model of total copy number is specified, it will automatically be converted to model of allele-specific copy number.

The copy number changes can be caused by at most five types of mutational events (segment duplication, segment deletion, chromosomal gain, chromosomal loss and whole genome doubling) along this tree for now. The mutation process is the same at all branches and across all sites.
The maximum copy number of the mutated genome is constrained by a user-specified parameter "--cn_max cn_max".


The mutation probability of each site is limited by its current copy number state.
Chromosomal gain is possible when the maximum copy number in the chromosome is smaller than the specified maximum copy number.
So when the specified mutation/duplication rate is high, the actual mutation rates along the lower branches gradually decrease due to copy number saturation.

The units of mutation rates are different at site (segment), chromosome, and whole genome levels. When computing the relative rates of a specific mutation type, they are summarized at different scales: total duplication/deletion rates obtained by summarizing over all sites, chromosome gain/loss rates obtained by summarizing over all chromosomes.

Please see run-sveta.sh to learn how to set different parameters

## Input
* --epop Ne: Ne is used to scale the tree height so that the branch length is measured in the unit of year, by multipling each branch length with Ne.

* --tiff delta_t: On the initial tree, the tip nodes have the same time. This parameter can be used to introduce different times at the tip nodes. The terminal branches are increased by random multiples of delta_t. The maximum multiple is the number of samples.

* --cn_max cn_max: The maximum copy number allowed in the program depends on the heap space.

* --cons 1/0: Whether or not the tree is constrained by patient age. If yes (1), the initial branch lengths will be adjusted by specified patient age so that the tree height is smaller than patient age.


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
The nodes in the tree are specified in a fixed order. The patient samples are indexed from 0 to n-1, where n is the number of samples. The normal sample has ID n. The root has ID n+1. For the internal nodes have ID from n+2 to 2n+1 from bottom to up.

File *-info.txt and *-mut.txt can be used to map mutations onto the tree.
The columns in *-mut.txt: sample_ID, edge_ID, muttype_ID, mut_btime, mut_etime, chr_haplotype, chr, seg_ID
For chromosome gain/loss, seg_ID is assigned to -1.
For whole genome doubling, chr is assigned to 0 and seg_ID is assigned to -1.



# Tree building with ML
<!-- ## How to prepare MP trees -->
## Input
The initial trees for tree searching can be obtained by maximum parsimony methods.
* (Required) A file containing copy numbers for all the samples, including the normal sample (*-cn.txt.gz or *-allele-cn.txt.gz)
* (Optional) A file containing the timing information of tip nodes (*-rel-times.txt)


There are 4 running modes in svtreeml.
* mode 0: building maximum likelihood tree from input copy numbers
* mode 1: a simple comprehensive test on a simulated tree
* mode 2: computing likelihood given a tree and its parameters (branch length and mutation rates)
* mode 3: computing maximum likelihood tree given a tree topolgy
The last three modes can be used to validate the computation of likelihood.

There are 3 tree searching method:
* exhaustive search (feasible for trees with fewer than 7 samples)
* hill climbing (only works for unconstrained case for now)
* genetic algorithm (may be slow, need improvement, not recommend to use)

Please see run-svtreeml.sh to learn how to set different parameters

There are four Markov models of evolution for building trees from the copy number profiles:
* model 0: Mk model (extension of JC69 model)
* model 1: bounded model of total copy number
* model 2: bounded model of allele-specific copy number
* model 3: independent Markov chain model (with 3 chains)

For tree recontruction with svtreeml on data with chromosome gain/loss and WGD, model 3 (model of independent Markov chains) should be used.
This model is DIFFERENT from the model used for simulating the data, which is usually model 2 (model of allele-specific copy number) in sveta.
Model 2 may also be used but there is a strong assumption of event order, that WGD is followed by chromosomal gain or loss and then segment duplication or deletion.


There are four important parameters for independent Markov chain model (model 3) to specify the number of states for each chain.
* max_wgd: the maximum number of WGD events in a sample. When it is 0,
* max_chr_change: the maximum number of chromosome gain/loss events accross all chromosomes in a sample. When it is 0, chromosome-level CNAs are not considered.
* max_site_change: the maximum number of segment duplication/deletion events accross all sites in a sample. When it is 0, segment-level CNAs are not considered.
* m_max: The maximum number of copies of a segment before chromosome gain/loss events.

Here, max_wgd, max_chr_change, and max_site_change determines the size of transition matrix of the 3 Markov chains for different levels of CNAs.
If one type of event is known to be missing in the data, please specify corresponding max_TYPE to be 0.

It is IMPORTANT to specify the appropriate values so that the probability of copy number state transition can be properly computed. Lower values will lead to wrong computation of tree likelihood. Higher values are perferred, but it will slow down the computation.

Since at most one WGD is allowed for each sample, max_wgd should be 1.
The maximum chromosome change (max_chr_change) should be 1 in most cases, but it should be increased accordingly if some chromosomes undergo more than one gain/loss events.
The segment-level CNAs are very likely to overlap, so max_site_change are often larger than 1.
Suppose there is no WGD and chromosome gain/loss, if the maximum copy number in the sample is 5, then max_site_change should be 3 (5 - 2).
You need adjust the values of max_site_change according to the input data.

You may check the copy number counts in the input data using similar command as below:
`less sim-data-1-cn.txt.gz | cut -f4 | sort | uniq -c`.



When estimating branch length, time constraints can be considered by specifying the following parameters.
* cons: When cons = 1, optimization is done with time constraints.


When building the tree, mutation rates can be estimated by specifying the following parameters.
* maxj: Whether or not to estimate mutation rate. When maxj = 1, mutation rates will be estimated. This is only reliable when the sampling times of tip nodes provide sufficient information.
* only_seg: When only_seg=1, only estimate segment duplication/deletion rates. Otherwise, the rates for chromosome gain/loss and WGD will be estimated.


## Output
* *-tree.txt: The reconstructed tree in tab-delimited format
* *-tree.nex: The reconstructed tree in NEWICK format, with branch length reprenting calender time





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

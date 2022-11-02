Copy Number Evolutionary Tree Analysis
=============

# Introduction
This repository contains a set of programs to simulate and build phylogenetic trees from copy number alterations (CNAs) in tumour genomes from multiple samples of a patient.

There are mainly 3 programs:
* cnets: simulating CNAs along a phylogenetic (coalescence) tree
* cnetml: building phylogenetic trees from copy numbers with maximum likelihood approach
* cnetmcmc (still under development): building phylogenetic trees from copy numbers with Bayesian MCMC approach


# Installation
This package is mostly written in C++. There are a few scripts written in R and Python, for plotting and text processing.

## Dependencies

* Required C/C++ libraries for compiling the main programs
  * [CMake](https://cmake.org/install/) is required for BFGS optimization
  * C libaries: [gsl](https://www.gnu.org/software/gsl/), [boost](https://www.boost.org/) (version >= 1.42)

* Required R libraries for postprocessing
  * `optparse`
  * `copynumber`
  * `tools`
  * `tidyverse`
  * `ggtree`
  * `ape`
  * `phangorn`


* Required Python libraries for postprocessing
  * newick2elist.py: [networkx](https://networkx.org)

### How to install CMake
#### Installing with package manager
On Mac:
```
brew install cmake
```

On Ubunbu:
```
sudo apt-get install cmake
```

#### Compiling from source
Get the latest “Unix/Linux Source” *.tar.gz file.
```
HOME=~/local
mkdir -p $HOME
tar -xf cmake*.tar.gz
cd cmake*
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

install.packages("optparse")
install.packages("tidyverse")
install.packages("ape")
install.packages("phangorn")


# Check installed packages
installed.packages()[, c("Package", "LibPath")]

```

## Building C++ source files
OpenMP is used to accelerate tree search in cnetml.
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
> bash run-cnets.sh

# Build a tree from copy number profile with maximum likelihood method
> bash run-cnetml.sh

# Build a tree from copy number profile with MCMC method
> bash run-cnetmcmc.sh
```

The most recent Mac switches to zsh. In that case, please replace `bash` with `zsh` in the commands above.



# Simulating data with cnets
cnets simulates structural variations that alter copy numbers along a phylogenetic tree of multiple samples.
Each tip in the tree corresponds to the genome of a sample, which is represented by a consecutive set of pre-specified sites.
Each site is considered as a segment of variable size or a bin of fixed size.
These sites can be seen as segments or bins obtained by segmentation methods or fixed-size sliding windows when calling copy numbers from real data.

The program generates a set of files, which record the simulated haplotype-specific and total copy number profiles, tree topology, mutations, and/or tip timing along the branches respectively.
The simulated copy number profiles and/or tip timing information served as input for the tree building methods.


The procedure of simulations is as follows:
1. Generate a random coalescence tree (with exponential growth) which represents the genealogy of tumour samples. Available trees can also be given as input, such as trees generated by [CoalEvol](https://code.google.com/archive/p/coalevol/). Add a branch of length zero which leads to a new top node representing a normal genome. The tree can be constrained by a specified patient age and the tips can be adjusted to reflect different sampling times。
		1) (optional, when tdiff > 0) Assign random times (in year) to the tips to create different sampling times.
		2) (optional, when constrained = 1) Rescale the tree by the specified patient age at first sampling time so that the tree height is no larger than the patient age at last sampling time.
2. Simulate copy number profiles on the tree according to a model of evolution.
3. Output result files

The maximum copy number of the mutated genome is constrained by a user-specified parameter "--cn_max cn_max".

There are three models of evolution for the copy number profiles:
<!-- * model 0: Mk model (extension of JC69 model) -->
* model 1: bounded model of total copy number
* model 2: bounded model of haplotype-specific copy number (default)
* model 3: infinite sites model

Mk model is an extension of JC69 model on total copy numbers, where duplication rate equals to deletion rate.

In bounded model, duplication rate may be different from deletion rate. Once a copy number becomes 0, it cannot be changed. Once a copy number becomes maximum (cn_max), it cannot be increased.

There are two ways of simulating mutations along a tree:
1. Simulating waiting times along a branch (default). A random waiting time is generated from the exponential distribution with mean `1/r`. Here, `r` is the total mutation rate across the genome, obtained by adding up the duplication and deletion rates across all sites in the genome, chromosomal gain or loss rates across all chromosomes and whole genome doubling rate. When a mutation is generated, its type is randomly chosen based on the relative rates of different types of mutational events.
2. Simulating sequences at the end of a branch. This is appropriate when the mutational events are not of interest. Only duplication and deletion are supported.

Note that when simulating waiting times, only haplotype-specific model is allowed. If model of total copy number is specified, it will automatically be converted to model of haplotype-specific copy number.
The copy number changes can be caused by at most five types of mutational events (duplication, deletion, chromosomal gain, chromosomal loss and whole genome doubling) along this tree for now. The mutation process is the same at all branches and across all sites.
The mutation probability of each site is limited by its current copy number state.
Chromosomal gain is ONLY possible when the maximum copy number in the chromosome is smaller than the specified maximum copy number.
The units of mutation rates are different at site, chromosome, and whole genome levels. When computing the relative rates of a specific mutation type, they are summarized at different scales: total duplication/deletion rates obtained by summarizing over all sites, chromosome gain/loss rates obtained by summarizing over all chromosomes.

Please see run-cnets.sh to learn how to set different parameters.

## Input
* `--epop Ne`: Ne (effective population size) is used to scale the tree height so that the branch length is measured in the unit of year, by multiplying each branch length with Ne.

* `--tdiff delta_t`: The tip nodes of an initial tree may have the same time. This parameter can be used to introduce different times at the tip nodes. The terminal branches are increased by random multiples of delta_t. The maximum multiple is the number of samples.

* `--cn_max cn_max`: The maximum copy number allowed in the program depends on the heap space.

* `--cons 1/0`: Whether or not the tree is constrained by patient age. If yes (1), the initial branch lengths will be adjusted by specified patient age so that the tree height is smaller than patient age.


## Output
* *-cn.txt.gz: The total copy number at each site for each sample, including the normal sample
* *-rel-times.txt: The sampling time of tip nodes
* *-allele-cn.txt.gz: The haplotype-specific copy number at each site for each sample
* *-rcn.txt.gz: The relative total copy number at each site for each sample
* *-allele-rcn.txt.gz: The processed haplotype-specific copy number by normalizing with a baseline at each site for each sample
* *-inode-cn.txt.gz: The copy number at each site for each internal node in the tree
* *-info.txt: The time of each node and the total number of mutations simulated on each branch of the tree, grouped by lineages of tip nodes.
* *-mut.txt: The list of simulated mutations on each branch of the tree.
* *-tree.txt: The simulated tree in tab-delimited format
* *-tree.nex: The simulated tree in NEWICK format, with branch length representing calendar time
* *-tree-nmut.nex: The simulated tree in NEWICK format, with branch length representing expected number of CNAs per site

File *-cn.txt.gz (*-rcn.txt.gz) can serve as the input to a tree building program that used absolute (relative) integer total copy number.

File *-allele-cn.txt.gz (*-allele-rcn.txt.gz) can serve as the input to a tree building program that used absolute (relative) integer haplotype-specific copy number.

File *-rel-times.txt can provide the timing information of tip nodes to allow estimation of divergence time and mutation rates.

Files *-tree.* provide the real tree, which can be used for measuring the accuracy of tree building programs.
The nodes in the tree are specified in a fixed order. The patient samples are indexed from 0 to n-1, where n is the number of samples. The normal sample has ID n. The root has ID n+1. For the internal nodes have ID from n+2 to 2n+1 from bottom to up.

File *-info.txt and *-mut.txt can be used to map mutations onto the tree.
The columns in *-mut.txt: sample_ID, edge_ID, muttype_ID, mut_btime, mut_etime, chr_haplotype, chr, seg_ID
For chromosome gain/loss, seg_ID is assigned to -1.
For whole genome doubling, chr is assigned to 0 and seg_ID is assigned to -1.



# Tree building with cnetml

cnetml is a new maximum likelihood method based on a novel evolutionary model of CNAs to infer phylogenies from spatio-temporal samples taken within a single patient. 


There are 4 running modes in cnetml.
* mode 0 (default): building maximum likelihood tree from input copy numbers, using "-b 1" for bootstrapping
* mode 1: a simple comprehensive test on a simulated tree
* mode 2: computing likelihood given a tree and its parameters (branch length and mutation rates)
* mode 3: computing maximum likelihood tree given a tree topology, using "-b 1" for bootstrapping
* mode 4: inferring marginal and joint ancestral states of a given tree from copy number profile

The last three modes can be used to validate the computation of likelihood.

Please see run-cnetml.sh to learn how to set different parameters.

When estimating branch length, time constraints (in longitudinal samples) can be considered by specifying the following parameters.
* `cons`: When cons = 1, optimization is done with time constraints on patient age and/or tip timing. Mutation rates can be estimated when there are considerate time differences among tips.

When building the tree, mutation rates can be estimated by specifying the following parameters.
* `estmu`: Whether or not to estimate mutation rate. When estmu = 1, mutation rates will be estimated. This is only reliable when the sampling times of tip nodes provide sufficient information (cons = 1 and dt > 0).
* `only_seg`: When only_seg = 1, only estimate duplication/deletion rates. Otherwise, the rates for chromosome gain/loss and WGD will be estimated.

Mutation rates are implicit parameters in the computing tree likelihood, so the reconstructed tree should be more accurate when cons = 1 and estmu = 1 given longitudinal samples, unless mutation rates are known as in simulated data.

There are 3 tree searching method:
* exhaustive search (feasible for trees with fewer than 7 samples)
* heuristic search (applicable for trees with at least 5 samples)
* genetic algorithm (may be slow, need improvement, deprecated)


There are four models of evolution for building trees from the copy number profiles:
* model 0: Mk model (deprecated)
* model 1: bounded model of total copy number (deprecated)
* model 2: bounded model of haplotype-specific copy number
* model 3: independent Markov chain model (with 3 chains, in development)

The first three models are the same as those for simulation.
The last one (model of independent Markov chains) should be used for tree reconstruction on data with chromosome gain/loss and WGD.
This model is DIFFERENT from the model used for simulating the data, which is usually model 2 (model of haplotype-specific copy number) in cnets.
Model 2 may also be used but there is a strong assumption of event order, that WGD is followed by chromosomal gain or loss and then duplication or deletion.

There are four important parameters for independent Markov chain model (model 3) to specify the number of states for each chain.
* `max_wgd`: the maximum number of WGD events in a sample. When it is 0,
* `max_chr_change`: the maximum number of chromosome gain/loss events across all chromosomes in a sample. When it is 0, chromosome-level CNAs are not considered.
* `max_site_change`: the maximum number of duplication/deletion events across all sites in a sample. When it is 0, site-level CNAs are not considered.
* `m_max`: The maximum number of copies of a segment before chromosome gain/loss events.

Here, max_wgd, max_chr_change, and max_site_change determine the size of transition matrix of the 3 Markov chains for different levels of CNAs.
If one type of event is known to be missing in the data, please specify corresponding max_TYPE to be 0.

It is IMPORTANT to specify the appropriate values so that the probability of copy number state transition can be properly computed. Lower values will lead to wrong computation of tree likelihood. Higher values are preferred, but it will slow down the computation.

Since at most one WGD is allowed for each sample, max_wgd should be 1.
The maximum chromosome change (max_chr_change) should be 1 in most cases, but it should be increased accordingly if some chromosomes undergo more than one gain/loss events.
The site-level CNAs are very likely to overlap, so max_site_change are often larger than 1.
Suppose there is no WGD and chromosome gain/loss, if the maximum copy number in the sample is 5, then max_site_change should be 3 (5 - 2).
You need to adjust the values of max_site_change according to the input data.
This can be done via the provided script with the following command:
`Rscript util/check_site_pattern.R -c sim1-cn.txt -t sim1-patterns.txt`

You may check the copy number counts in the input data using similar command as below:
`less sim-data-1-cn.txt.gz | cut -f4 | sort | uniq -c`.


<!-- ## How to prepare MP trees -->

## Input
* (Required) A file containing integer absolute/relative copy numbers for all the patient samples and/or the normal sample (*-cn.txt.gz or *-allele-cn.txt.gz). When the input copy numbers are relative with normal copy being 0 as those output by [CGHcall](https://bioconductor.org/packages/release/bioc/html/CGHcall.html), please specifiy it with option "--is_rcn 1". When the input copy numbers are haplotype-specific which have been scaled relative to ploidy or not, please specifiy it with option "--is_total 0 --is_rcn 0".

   Either compressed file or uncompressed file is fine.
   There need to be at least four columns, separated by space, in this file: sample_ID, chr_ID, site_ID, CN.
   Note that there should be no header names in this file.
   Each column is an integer.
   The sample_ID has to be __ordered__ from 1 to the number of patient samples.
   The chr_ID and site_ID together determine a unique site along the genome of a sample, __ordering__ from 1 to the largest number.
   The site_ID can be consecutive numbers from 1 to the total number of sites along the genome, or consecutive numbers from 1 to the total number of sites along each chromosome of the genome.
   For haplotype-specific CN, there need to be at least five columns, with the last two being cnA, cnB.
   If the total CN is larger than the specified maximum CN allowed by the program,
   the total CN will be automatically decreased to the maximum CN when the input is total CN and the program will exit when the input is haplotype-specific CN.

* (Optional) A file containing the timing information of patient samples (*-rel-times.txt).

  There need to be three columns, separated by space, in this file:
  sample_ID, time relative to 1st sample in year (float number), patient age at the time of sampling (integer).
  The sample_ID has to be ordered from 1 to n (the number of patient samples).

* (Optional) A folder containing initial trees for tree searching, which may obtained by parsimony-based methods.

    It is recommended when the tree is large to avoid local optima.
  
## Output
* *-tree.txt: The reconstructed tree in tab-delimited format.
* *-tree.nex: The reconstructed tree in NEWICK format, with branch length representing calendar time.
* *-tree.nmut.nex: The reconstructed tree in NEWICK format, with branch length representing number of mutations.
* *-segs.txt: The file with the postprocessed copy number matrix for tree building, with columns being chr_D, start_bin_ID, end_bin_ID, start_segment_ID, end_segment_ID, and the copy numbers for each sample. It will be generated when reading the input copy number file.

## Note 
Please ensure the input file exist and their names are correct, or else there may be an error of "Segmentation fault (core dumped)"!


## Ancestral state reconstruction
Only the implementation under the bounded model of haplotype-specific copy number has been fully tested so far.

## Input
* (Required) A file containing integer absolute copy numbers for all the patient samples and/or the normal sample (*-cn.txt.gz or *-allele-cn.txt.gz).
* (Required) *-tree.txt: A tree file in tab-delimited format.
* (Required) --dup_rate dup_rate --del_rate del_rate: The mutation rates of the input tree.
* (Optional) A file containing the timing information of patient samples (*-rel-times.txt).

## Output
* *.(mrca|joint).cn: The reconstructed copy numbers in tab-delimited format for the most recent common ancestor node of all tumour samples (mrca) and all internal nodes (joint) respectively. When the input are total copy numbers, there are four columns:  node ID, chromosome ID, site ID, CN. When the input are haplotype-specific copy numbers, there are five columns: node ID, chromosome ID, site ID, cnA, cnB.
* *.mrca.state: A tab-delimited file containing the posterior probability of each possible copy number state on a unique variant site (containing at least one atypical copy number across all samples) for MRCA node. The columns are: node ID, site ID (chromosomeID_siteID), probability_stateID. The copy number state ID is the same as the state index in the rate matrix of the Markov model.
* *.joint.state: A tab-delimited file containing the possible copy number state on a unique variant site for all internal nodes. The columns are: node ID, site ID (chromosomeID_siteID), cn_stateID.



# Tree building with cnetmcmc

Only basic MCMC algorithm is implemented here and not comprehensively tested.

Please see run-cnetmcmc.sh to learn how to set different parameters.

There are two running modes depending on whether a reference tree is provided or not.
With a reference tree, the tree topolgy is fixed.

## Input
* (Required) A file containing copy numbers for all the samples, including the normal sample (*-cn.txt.gz or *-allele-cn.txt.gz)
* (Optional) A file containing the timing information of tip nodes (*-rel-times.txt)
* (Optional) A configuration file which sets most input parameters (mcmc.cfg)


## Output
There are two output files in a format similar to that of MrBayes:
* *.p, which records the traces of parameters
* *.t, which records the sampled trees

<!-- ## How to analyze the results of cnetmcmc -->
The output can be analyzed by [RWTY](https://github.com/danlwarren/RWTY). Please see script util/check_convergence.R for reference.

*.p can be imported into [Tracer](https://beast.community/tracer) to check the convergence of the chains.

*.t can be analyzed by [TreeAnnotator](https://beast.community/treeannotator) to get a summary tree (maximum credibility tree).

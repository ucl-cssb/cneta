#!/usr/bin/bash

# This script is used to run program sveta, which can generate random colesence trees and mutations along the branches.

seed=1218794168  # used for reproducing the results
# Whether or not to print debug information
verbose=0

# The number of regions to sample. Output will be Nr+1 including germline
Ns=3

# The number of simulations. Number of tumours to simulate
Nsim=1

# Maximum copy number allowed (failed when > 13)
cn_max=4

mode=1  # 0: Simuting genome in fix-sized bins, 1: Simulating genome in segments of random size
seg_max=1000 # Maximum number of segments on the genome
fix_seg=1 # Whether or not to fix the number of segments to be seg_max

model=2 # Model of evolution. 1: bounded model, 0: JC69 model, 2: allele-specific model

method=0  # Simulating method. 0: waiting time, 1: direct sequences

# Whether the tree is constrained by age or not
cons=0
# Age of the patient
age=60
# time step for simulating different sampling times
dt=0


# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.01
r2=0.01
r3=0
r4=0
r5=0
# mean of exponential distributions of duplication and deletion size in bins
s1=5
s2=5


# The output directory
# dir="./example/method${method}"
dir="./example"
if [[ ! -d $dir ]]; then
  mkdir -p $dir
fi

tree_file=""  # The input tree file. If given, mutations will be generated along this tree
# effective population size for the coalescence tree
Ne=9000000
beta=1.563e-3  # exponential growth rate
gtime=0.002739726 # generation time in year

# Prefix of output file. Please set it to be "" when generating multiple samples to avoid overwriting previous results
# prefix=sim-data-1
prefix=""


code/sveta --tree_file "$tree_file" -o $dir -r $Ns -n $Nsim --mode $mode --method $method --fix_seg $fix_seg --seg_max $seg_max --cn_max $cn_max --dup_rate $r1 --del_rate $r2 --chr_gain $r3 --chr_loss $r4 --wgd $r5 --dup_size $s1 --del_size $s2 -e $Ne -b $beta -t $dt --verbose $verbose --constrained $cons --model $model -p "$prefix" --age $age --seed $seed > $dir/std_sveta_cons"$cons"_model"$model"_method"$method"
#

# Plot all simulated trees
Rscript ana/plot-trees-all.R -d $dir -b 0 -t "all" -l "xlim"  # >& /dev/null
# Plot simulated tree with the number of mutations on the branch
Rscript ana/plot-trees-all.R -d $dir -b 1 -t "all" # >& /dev/null
Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata # >& /dev/null

# Plot a single tree
# prefix=sim-data-1
# # Simulted sampling time informaton
# tfile=$dir/${prefix}-rel-times.txt
# ofile=$dir/"$prefix"-tree.txt
# cfile=$dir/"$prefix"-cn.txt.gz
# Rscript ana/plot-trees-all.R -f $ofile -b 0 -t "single" -l "xlim" --time_file $tfile # >& /dev/null
# # Plot simulated tree with the number of mutations on the branch
# Rscript ana/plot-trees-all.R -f $ofile -b 1 -t "single"  # >& /dev/null
# Rscript ana/plot-cns.R -f $cfile -b ana/bin_locations_4401.Rdata # >& /dev/null

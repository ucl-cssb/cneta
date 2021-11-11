#!/usr/bin/env bash

# This script is used to run program sveta, which can generate a random colesence tree (of tumor samples from a single patient) and structual variations along the tree branches.

seed=$RANDOM  # used for reproducing the results
verbose=0   # Whether or not to print debug information
Nsim=1  # The number of simulations. Number of patients to simulate


####################### Parameters related to tree generation ##################
Ns=3 # The number of tumor regions to sample. Output will be Ns+1 including germline normal region
tree_file=""  # The input tree file. If given, mutations will be generated along this tree
Ne=9000000  # effective population size for the coalescence tree
beta=1.563e-3  # exponential growth rate
gtime=0.002739726 # generation time in year
age=60  # Age of the patient at first sample
cons=1  # Whether the tree height is constrained by age or not
dt=2  # time step for simulating different sampling times
################################################################################


####################### Parameters related to mutation generation ##############
model=2 # Model of evolution. 0: Mk, 1: one-step bounded (total), 2: one-step bounded (allele-specific), 3: Poisson
mode=1  # 0: Simuting genome in fix-sized bins, 1: Simulating genome in segments of random size
seg_max=1000  # Maximum number of segments on the genome
fix_nseg=1  # Whether or not to fix the number of segments to be seg_max
method=1  # Simulating method. 0: waiting times, 1: sequences

cn_max=4  # Maximum copy number allowed (program failed due to memory issue when > 13)
# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.001
r2=0.001
r3=0
r4=0
r5=0
# mean of exponential distributions of duplication and deletion size in bins
s1=5
s2=5
################################################################################


####################### Set output directory and run simulation ################
# The output directory
# dir="./example/method${method}"
dir="./example"
if [[ ! -d $dir ]]; then
  mkdir -p $dir
fi

# Prefix of output file. Please set it to be "" when generating multiple samples to avoid overwriting previous results
# prefix=sim-data-1
prefix=""

echo "Start running sveta"

echo "seed $seed"  > $dir/std_sveta_cons${cons}_model${model}_method${method}
# valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v
/usr/bin/time ./code/sveta --tree_file "${tree_file}" -p "${prefix}" -o $dir -r $Ns -n $Nsim --mode $mode --method $method --fix_nseg $fix_nseg --seg_max $seg_max --cn_max $cn_max --dup_rate $r1 --del_rate $r2 --chr_gain $r3 --chr_loss $r4 --wgd $r5 --dup_size $s1 --del_size $s2 -e $Ne -b $beta --gtime $gtime -t $dt --verbose $verbose --constrained $cons --model $model --age $age --seed $seed >> $dir/std_sveta_cons${cons}_model${model}_method${method}

wait

echo "Finish running sveta"
################################################################################


####################### Plot simulated tree and copy numbers (optional) ########
# Plot all simulated trees
Rscript ana/plot-trees-all.R -d $dir -b 0 -t "all" -l "xlim"  # >& /dev/null
# Plot simulated tree with the number of mutations on the branch
Rscript ana/plot-trees-all.R -d $dir -b 1 -t "all" -l "xlim" # >& /dev/null
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
################################################################################

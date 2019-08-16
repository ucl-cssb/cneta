#!/usr/bin/bash

# seed=847596759  # used for reproducing the results
# Whether or not to print debug information
verbose=0

# The output directory
dir="./example/"
mkdir $dir
# Prefix of output file. Please set it to be "" when generating multiple samples to avoid overwriting previous results
prefix=sim-data-1
# prefix=""

# The number of regions to sample. Output will be Nr+1 including germline
Ns=4
# The number of simulations. Number of tumours to simulate
Nsim=1

# Whether the tree is constrained by age or not
cons=1
# Model of evolution. 1: bounded model, 0: JC69 model
model=1
# Age of the patient
age=60
# Maximum copy number allowed
cn_max=4
# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.00001
r2=0.00002
r3=0
r4=0
r5=0
# mean of exponential distributions of duplication and deletion size in bins
s1=50
s2=50
# effective population size for the coalescence tree
Ne=10
# time step for simulating different sampling times
dt=2


code/sveta -o $dir -r $Ns -n $Nsim --cn_max $cn_max --dup_rate $r1 --del_rate $r2 --chr_gain $r3 --chr_loss $r4 --wgd $r5 --dup_size $s1 --del_size $s2 -e $Ne -t $dt --verbose $verbose --constrained $cons --model $model -p "$prefix" --age $age > $dir/std_sveta_cons"$cons"_model"$model"
# --seed $seed
# Simulted sampling time informaton
tfile=$dir/${prefix}-rel-times.txt

# Plot all simulated trees
# Rscript ana/plot-trees-all.R -d $dir -b 0 -t "all" -l "xlim"  # >& /dev/null
# # Plot simulated tree with the number of mutations on the branch
# Rscript ana/plot-trees-all.R -d $dir -b 1 -t "all" -l "xlim" # >& /dev/null
# Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata # >& /dev/null

# Plot a single tree
ofile=$dir/"$prefix"-tree.txt
Rscript ana/plot-trees-all.R -f $ofile -b 0 -t "single" -l "age" --time_file $tfile # >& /dev/null
# Plot simulated tree with the number of mutations on the branch
Rscript ana/plot-trees-all.R -f $ofile -b 1 -t "single"  # >& /dev/null
Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata -p "${prefix}-cn.txt.gz" # >& /dev/null

#!/usr/bin/bash

# output directory
dir="./test/"
mkdir $dir

# seed=165873033

# number of regions. Output will be Nr+1 including germline
Ns=5
# number of samples. Number of tumours to simulate
Nsim=1
# Whether the tree is constrained by age or not
cons=1
# Model of evolution. 1: bounded model, 0: JC69 model
model=1
# Age of the patient
age=50

# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.00004
r2=0.00002
r3=0.00001
r4=0.00001
r5=0.000001

# mean of exponential distributions of duplication and deletion size in bins
s1=20
s2=20

# effective population size
Ne=1.0
# observed time step
dt=0.5

# Whether or not to print debug information
verbose=0

# Prefix of output file
prefix=""
# prefix=sim-data-N${Ns}-cons${cons}-model${model}-1


code/sveta -o $dir -r $Ns -n $Nsim --dup_rate $r1 --del_rate $r2 --chr_gain $r3 --chr_loss $r4 --wgd $r5 --dup_size $s1 --del_size $s2 -e $Ne -t $dt --verbose $verbose --constrained $cons --model $model -p "$prefix" --age $age > $dir/std_test_sveta_cons"$cons"_model"$model"
# --seed $seed

# Plot all simulated trees
Rscript ana/plot-trees-all.R -d $dir -b 0 -t "all" # >& /dev/null
# Plot simulated tree with the number of mutations on the branch
Rscript ana/plot-trees-all.R -d $dir -b 1 -t "all" # >& /dev/null
Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata # >& /dev/null

# Plot a single tree
# ofile=$dir/"$prefix"-tree.txt
# Rscript ana/plot-trees-all.R -f $ofile -b 0 -t "single" -l "xlim" # >& /dev/null
# # Plot simulated tree with the number of mutations on the branch
# Rscript ana/plot-trees-all.R -f $ofile -b 1 -t "single"  # >& /dev/null
# Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata -p "${prefix}-cn.txt.gz" # >& /dev/null

#!/usr/bin/bash

# Use MCMC to infer tree
# With configuration file

is_total=1 # If yes, the input is total copy number

# input parameters
dir="./example"   # The output directory
# prefix=sim-data-N5-cons1-model1-1
prefix=sim-data-1
input=$dir/"$prefix"-cn.txt.gz
# input=$dir/"$prefix"-allele-cn.txt.gz
times=$dir/"$prefix"-rel-times.txt
# rtree=$dir/"$prefix"-tree.txt

rtreefile=""    # Optional reference tree
nsample=3

odir=$dir/mcmc
mkdir -p $odir

suffix=$prefix
config_file=/mnt/d/Gdrive/git/sveta/mcmc.cfg

# Run multiple chains
for i in {1..2..1}
do
  # seed=$i
  echo $i
  trace_param_file=$odir/trace-mcmc-params_"$suffix"_${i}.txt
  trace_tree_file=$odir/trace-mcmc-trees_"$suffix"_${i}.txt
  sum_tree_file=$odir/sum-mcmc-trees_"$suffix"_${i}.txt

  code/svtreemcmc --nsample $nsample -c $input -t $times --rtree "$rtree" --trace_param_file $trace_param_file --trace_tree_file $trace_tree_file --config_file $config_file > $odir/std_mcmc_"$suffix"_${i}

  # Summarize the sampled trees into a maximum credibility tree with median heights
  # treeannotator -burnin 10 -heights median $trace_tree_file $sum_tree_file
done

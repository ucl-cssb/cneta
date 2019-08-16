#!/usr/bin/bash

# Use MCMC to infer tree
# With configuration file

# input parameters
dir=./example
# prefix=sim-data-N5-cons1-model1-1
prefix=sim-data-1
input=$dir/"$prefix"-cn.txt.gz
times=$dir/"$prefix"-rel-times.txt
# rtree=$dir/"$prefix"-tree.txt
rtreefile=""
nsample=4

odir=$dir/mcmc
mkdir -p $odir

suffix=$prefix
config_file=/mnt/d/Gdrive/git/sveta/mcmc.cfg


for i in {1..1..1}
do
  # seed=$i
  echo $i
  trace_param_file=$odir/trace-mcmc-params_"$suffix"_${i}.txt
  trace_tree_file=$odir/trace-mcmc-trees_"$suffix"_${i}.txt
  code/svtreemcmc --nsample $nsample -c $input -t $times --rtree "$rtree" --trace_param_file $trace_param_file --trace_tree_file $trace_tree_file --config_file $config_file > $odir/std_mcmc_"$suffix"_${i}
done

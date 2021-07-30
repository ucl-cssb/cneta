# This script is used to run program svtreemcmc.

# Use MCMC to infer tree
# With configuration file
config_file=./mcmc.cfg


####################### Parameters related to input  ###########################
Ns=3  # The number of regions

# input parameters
dir="./example"   # The output directory
# prefix=sim-data-N5-cons1-model1-1
prefix=sim-data-1
is_total=1 # If yes, the input is total copy number
input=$dir/"$prefix"-cn.txt.gz
# input=$dir/"$prefix"-allele-cn.txt.gz

# The input file of sample timings (optional, required for estimating mutation rates)
times=$dir/"$prefix"-rel-times.txt
if [[ ! -f $times ]]; then
  times=""
fi

# Optional reference real tree
# rtree=$dir/"$prefix"-tree.txt
rtree=""

init_tree=0     # Types of starting tree for MCMC. 0: random tree, 1: provided tree, 2: random tree with the same topology as real tree
file_itree=$dir/"$prefix"-tree.txt
################################################################################


####################### Set output directory and run the program ###############
# output parameters
odir=$dir/mcmc
if [[ ! -d $odir ]]; then
  mkdir -p $odir
fi
suffix=$prefix

echo "Start running svtreemcmc"

# Run multiple chains
for i in {1..2}
do
  # seed=$i
  echo "chain $i"
  # Use the suffix of MrBayes for compatibility with RWTY
  trace_param_file=$odir/mcmc_"$suffix"_${i}.p
  trace_tree_file=$odir/mcmc_"$suffix"_${i}.t
  sum_tree_file=$odir/sum-mcmc_"$suffix"_${i}.txt

  code/svtreemcmc -s $Ns --is_total $is_total -c $input -t "$times" --rtree "$rtree" --trace_param_file $trace_param_file --trace_tree_file $trace_tree_file --config_file $config_file --init_tree $init_tree --file_itree $file_itree > $odir/std_mcmc_"$suffix"_${i}

  # Summarize the sampled trees into a maximum credibility tree with median heights
  # treeannotator -burnin 10 -heights median $trace_tree_file $sum_tree_file
done

echo "Finish running svtreemcmc"

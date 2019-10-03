#!/usr/bin/bash

seed=349723784 # Setting seed for reproductive results
verbose=1  # Whether or not to print debug information

idir="./example/method0"   # The output directory
prefix=sim-data-1   # The prefix of output files
dir=$idir
# dir=$idir/total
mkdir -p $dir

# The input parametersW
# input=$idir/"$prefix"-cn.txt.gz
input=$idir/"$prefix"-allele-cn.txt.gz
times=$idir/"$prefix"-rel-times.txt
init_tree=0  # 0: Random coalescence tree, 1: Maximum parsimony tree
dir_itrees=$idir/itrees  # The directory containing the initial trees (in TXT format)
tree_file=$idir/"$prefix"-tree.txt   # Only required when the tree is to be given as input

# effective population size for random initial coalescence tree
Ne=90000
beta=1.563e-3  # exponential growth rate
gtime=0.002739726 # generation time in year

# Parameters that should be consistent with the simulation or real data
Ns=3  # The number of regions
cn_max=4  # Maximum copy number allowed
model=2  # Model of evolution.  0: Mk model, 1: bounded model, 2: allele-specific model
age=60  # Age of the patient
# Specify the (initial) mutation rates when they are fixed
r1=0.001
r2=0.001
r3=0
r4=0
r5=0

is_total=0
is_bin=0
incl_all=0
mode=0  # Running mode. 0: build maximum likelihood tree; 1: test; 2: compute likelihood; 3: compute maximum likelihood
cons=0  # Whether or not the tree is constrained by patient age
maxj=0  # Whether or not to estimate mutation rate
only_seg=1  # 1: only consider segment duplication/deletion; 0: otherwise
mu=0  # Used in model 0

# Parameters for optimization and tree search algorithms
correct_bias=1  # Whether or not to correct acquisition bias in likelihood computation
opt=1   # 1: L-BFGS-B method; 0: simplex method
tolerance=1e-2
tree_search=2 # 0: genetic algorithm, 1: hill climbing, 2: exhaustive search
Npop=100
Ngen=20
Nstop=5  #


if [[ $mode -eq 0 ]]; then
  suffix=o"$opt"-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  code/svtreeml -c $input -t $times --tree_file $tree_file --is_total $is_total --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj  --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode --seed $seed > $dir/std_svtreeml_"$suffix"
  #

  Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim" --time_file $times  #>& /dev/null
  Rscript ana/plot-trees-all.R -f $mltree -b 1 -t "single" -l "xlim" --time_file $times  #>& /dev/null

  bootstap=0
  if [[ $bootstap -eq 1 ]]; then
    # Do bootstapping
    bdir=$dir/bootstrap
    mkdir -p $bdir
    for i in {1..10..1}
    do
      echo $i
      ofile=$bdir/MaxL-"$suffix"-btree-$i.txt

      code/svtreeml -c $input -t $times --tree_file $tree_file --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $ofile -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj  --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2  --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode -b 1 > $bdir/std_svtreeml_"$suffix"-btree-$i
    done
    # Draw the ML tree with bootstapping support
    Rscript ana/plot-trees-all.R -s $bdir -f $mltree -o $dir/MaxL-tree-"$suffix"-bootstrap.pdf -t "bootstrap" -l "age" --time_file $times -p "MaxL-"$suffix"-btree-*txt"
  fi

elif [[ $mode -eq 1 ]]; then
  suffix=sim1-m$model-test-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  code/svtreeml -c $input -t $times --is_bin $is_bin --incl_all $incl_all --tree_file $tree_file -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose > $dir/std_svtreeml_"$suffix"

  Rscript ana/plot-trees-all.R -d $dir/ -b 0 -t "all" -l "xlim" --time_file $times -p "sim-data-*-tree*_bfgs.txt"
  Rscript ana/plot-trees-all.R -d $dir/ -b 1 -t "all" -l "xlim" --time_file $times -p "sim-data-*-tree*_bfgs.txt"

elif [[ $mode -eq 2 ]]; then
  # In this mode, the mutation rates have to be specified
  suffix=m$model-"$cons""$maxj"-mode"$mode"-"$prefix"-all"$1"

  code/svtreeml -c $input -t $times --is_bin $is_bin --incl_all $incl_all -s $Ns --tree_file $tree_file -d $model --cn_max $cn_max --only_seg $only_seg --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_svtreeml_"$suffix"

else
  # In this mode, the mutation rates can be specified or not
  suffix=m$model-o"$opt"-"$cons""$maxj"-mode"$mode"-"$prefix"-all"$1"
  mltree=$dir/MaxL-"$prefix"-"$suffix".txt

  code/svtreeml -c $input -t $times --is_bin $is_bin --incl_all $incl_all -s $Ns --tree_file $tree_file -o $mltree -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_svtreeml_"$suffix"

  Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim" --time_file $times  #>& /dev/null
fi


# Estimating mutation rates given the tree
# mutree=$dir/results-mu-tree-"$suffix".txt
# code/svtreemu -c $input -t $times -p  $mltree -s $Ns -x $mut -l 0.5  -d $model -o $mutree -m 2000 -r 0.01 > $dir/std_svtreemu_"$suffix"
# #Rscript ana/plot-trees-single.R $dir/MaxL-mu-tree.txt 0 #>& /dev/null
# Rscript ana/plot-trees-all.R -f $mutree -b 0 -t "single" #>& /dev/null

#!/usr/bin/env bash

# This script is used to run program cnetml.

seed=$RANDOM # Setting seed for reproductive results
verbose=1  # Whether or not to print debug information

# Running mode. 0: build maximum likelihood tree; 1: test; 2: compute likelihood; 3: compute maximum likelihood; 4: infer ancestral state
mode=0
plot=0

####################### Parameters related to input  ###########################
idir="./example"   # The input directory
prefix=sim-data-1   # The prefix of input files
# prefix=sim-data-"$1"

# The input parameters related to copy numbers
input=$idir/"$prefix"-cn.txt.gz
# input=$idir/"$prefix"-allele-cn.txt.gz
if [[ ! -f $input ]]; then
  echo "Input file $input does not exit!"
  exit 0
fi
seg_file=$idir/"$prefix"-segs.txt

is_total=1    # If yes (1), the input is total copy number; or else (0), the input is haplotype-specific
# input=$idir/"$prefix"-allele-cn.txt.gz
# whether or not the input copy number is for each bin. If not, the input copy number is read as it is. Or else, consecutive bins will be merged
is_bin=0
# whether or not to include all the input copy numbers without further propressing
incl_all=1
cn_max=4  # Maximum copy number allowed

# The input parameters that should be consistent with the simulation or real data
Ns=3 # The number of samples for a patient
cons=1  # Whether or not the tree is constrained by patient age

# The input file of sample timings (optional, required for estimating mutation rates)
times=$idir/"$prefix"-rel-times.txt
if [[ ! -f $times ]]; then
  times=""
fi

model=2  # Model of evolution.  0: Mk model, 1: bounded model for total copy number, 2: bounded model for haplotype-specific copy number, 3: independent Markov chain model
# Parameters for independent Markov chain model
max_wgd=0
max_chr_change=0
max_site_change=3
m_max=1

estmu=1  # Whether or not to estimate mutation rate
only_seg=1  # 1: only estimate site duplication/deletion rates; 0: otherwise

if [[ $times == "" ]]; then
  estmu=0
  cons=0
fi
################################################################################


######### Parameters for optimization and tree search algorithms  ##############
correct_bias=1  # Whether or not to correct acquisition bias in likelihood computation
opt=1   # 1: L-BFGS-B method; 0: simplex method
tolerance=1e-2
miter=2000

tree_search=2  # 0: genetic algorithm, 1: hill climbing, 2: exhaustive search
# Parameters used in genetic algorithm
Npop=100  # size of initial trees for hill climbing
Ngen=20
Nstop=5
# whether or not to do reduced NNI while doing hill climbing NNIs
speed_nni=1

# Parameters about starting tree during search
init_tree=0  # 0: Random coalescence tree, 1: Maximum parsimony tree
dir_itrees=$idir/itrees  # The directory containing the initial trees (in TXT format)
tree_file=$idir/"$prefix"-tree.txt   # Only required when the tree is to be given as input
# tree_file=""
# effective population size for random initial coalescence tree
Ne=90000
beta=1.563e-3  # exponential growth rate
gtime=0.002739726 # generation time in year


# Specify the (initial) mutation rates when they are fixed
r1=0.001
r2=0.001
r3=0
r4=0
r5=0
mu=0  # Used in model 0
################################################################################


####################### Set output directory and run the program ###############
echo "Start running cnetml"

dir=$idir
if [[ ! -d $dir ]]; then
  mkdir -p $dir
fi
bsdir1=$dir/bootstrap
bsdir2=$dir/bootstrap_fixt

if [[ $mode -eq 0 ]]; then
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  echo "seed $seed" > $dir/std_cnetml_"$suffix"
  # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v
  /usr/bin/time code/cnetml -c $input --seg_file $seg_file -t "$times" --tree_file "$tree_file" --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode --speed_nni $speed_nni --seed $seed 2>&1 >> $dir/std_cnetml_"$suffix"
  # #
  #
  echo "Finish running cnetml"

  if [[ $plot -eq 1 ]]; then
    Rscript util/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim" --time_file "$times"  #>& /dev/null
  fi
  # Evaluate the estimation error
  # cmp_plot=$dir/cmp_plot-"$suffix".pdf
  # cmp_dist=$dir/cmp_dist-"$suffix".txt
  # Rscript test/tree_comparison/compare_trees.R -r "$tree_file" -i $mltree -o $cmp_plot -s $cmp_dist -p 1

  bootstrap=0
  if [[ $bootstrap -eq 1 ]]; then
    echo "Run bootstrapping"
    # Do bootstrapping
    mkdir -p $bsdir1
    for i in {1..10..1}
    do
      echo $i
      ofile=$bsdir1/MaxL-"$suffix"-btree-$i.txt

      code/cnetml -c $input -t "$times" --tree_file "$tree_file" --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $ofile -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu  --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2  --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode -b 1 > $bsdir1/std_cnetml_"$suffix"-btree-$i
    done
    # Draw the ML tree with bootstrapping support
    if [[ $plot -eq 1 ]]; then
      Rscript util/plot-trees-all.R -s $bsdir1 -f $mltree -o $dir/MaxL-tree-"$suffix"-bootstrap.pdf -t "bootstrap" -l "age" --time_file "$times" -p "MaxL-"$suffix"-btree-*txt"
    fi
  fi

elif [[ $mode -eq 1 ]]; then
  suffix=sim1-m$model-test-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  code/cnetml -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all --tree_file "$tree_file"  --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose > $dir/std_cnetml_"$suffix"

elif [[ $mode -eq 2 ]]; then
  # In this mode, the mutation rates have to be specified
  suffix=m$model-"$cons""$estmu"-mode"$mode"-"$prefix"

  code/cnetml -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns --tree_file "$tree_file" --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change -d $model --cn_max $cn_max --only_seg $only_seg --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode  --seed $seed > $dir/std_cnetml_"$suffix"

elif [[ $mode -eq 3 ]]; then
  # In this mode, the mutation rates can be specified or not
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt
  mltree2=$dir/MaxL-"$suffix"-fixt.txt
  echo "Estimating the branch lengths (and mutation rates) of the given tree $mltree"
  if [[ ! -f $mltree ]]; then
  echo "$mltree does not exist!"
  exit
  fi

  # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v
  /usr/bin/time code/cnetml --tree_file "$mltree" -o $mltree2 -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_cnetml_"$suffix"-mode"$mode"

  # plot the new tree
  if [[ $plot -eq 1 ]]; then
    Rscript util/plot-trees-all.R -f $mltree2 -b 0 -t "single" -l "xlim" --time_file "$times"  #>& /dev/null
  fi

  bootstrap=1
  if [[ $bootstrap -eq 1 ]]; then
    echo "Run bootstrapping"
    # Do bootstrapping
    mkdir -p $bsdir2

    for i in {1..10..1}
    do
      echo $i
      ofile=$bsdir2/MaxL-"$suffix"-btree-$i.txt

      # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v /usr/bin/time
      code/cnetml --tree_file "$mltree" -o $ofile -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns  --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change  -d $model --cn_max $cn_max --only_seg $only_seg --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode -b 1 > $bsdir2/std_cnetml_"$suffix"-btree-$i
    done

    # Draw the ML tree with confidence intervals
    ci_dup=`Rscript util/compute_ci.R $fdup 4`
    ci_del=`Rscript util/compute_ci.R $fdel 4`

    echo $ci_dup
    echo $ci_del

    if [[ $plot -eq 1 ]]; then
      Rscript util/plot-trees-all.R -s $bsdir1 -f $mltree2 -o $dir/MaxL-tree-"$suffix"-fixt-bootstrap.pdf -t "bootstrap" -l "ci" --time_file "$times" -p "MaxL-"$suffix"-btree-*txt" --bstrap_dir2 $bsdir2 --title "duplication rate $ci_dup; deletion rate $ci_del"
    fi
  fi

elif [[ $mode -eq 4 ]]; then  # Inferring ancestral states of a given tree from copy number profile
  # In this mode, the mutation rates have to be specified
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt
  echo "Inferring ancestral states of the given tree $mltree"

  code/cnetml -c $input -t "$times" --tree_file "$mltree" -o "$mltree" --is_bin $is_bin --incl_all $incl_all -s $Ns --is_total $is_total --max_wgd $max_wgd --max_chr_change $max_chr_change --max_site_change $max_site_change -d $model --cn_max $cn_max --only_seg $only_seg --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode  --seed $seed > $dir/std_cnetml_state_"$suffix"
fi

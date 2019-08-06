#!/usr/bin/bash

# input parameters
# prefix=sim-data-N"$Ns"-cons"$cons"-model"$model"-1
prefix=sim-data-1
dir=./example
input=$dir/"$prefix"-cn.txt.gz
times=$dir/"$prefix"-rel-times.txt
# tree_file=$dir/"$prefix"-tree.txt
# tree_file=$dir/AllTreesNr5_"$1".txt
dir_itrees=$dir/itrees  # The directory containing the initial trees (in TXT format)

Ns=10  # number of regions

mode=1 # Running mode. 0: test; 1: build maximum likelihood tree; 2: compute likelihood; 3: compute maximum likelihood

model=1 # Model of evolution. 1: bounded model, 0: JC69 model
cn_max=4
opt=1 # 1: L-BFGS-B method; 0: simplex method
cons=1  # Whether or not the tree is constrained by patient age
maxj=1 # Whether or not to estimate mutation rate
correct_bias=1  # Whether or not to correct acquisition bias
tree_search=1 # 1: hill climbing; 0: genetic algorithm
init_tree=1  # 0: Random coalescence tree, 1: Maximum parsimony tree

# Age of the patient
age=60

# Whether or not to print debug information
verbose=0

# Specify the mutation rates when they are fixed
r1=0.00001
r2=0.00002
r3=0
r4=0
r5=0

mu=0  # Used in model 0 (JC69 model)

# seed=1215208885 # Setting seed for reproductive results

# tree search algorithm
Npop=100
Ngen=100
Nstop=10


if [[ $mode -eq 1 ]]; then
  suffix=o"$opt"-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o $mltree -d $model --cn_max $cn_max --tree_search $tree_search --init_tree $init_tree --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj  --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2  --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_svtreeml_"$suffix"
  # --seed $seed

  Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim"  #>& /dev/null

  bootstap=0
  if [[ $bootstap -eq 1 ]]; then
    # Do bootstapping
    bdir=$dir/bootstrap
    mkdir -p $bdir
    for i in {1..10..1}
    do
      echo $i
      ofile=$bdir/results-maxL-"$suffix"-btree-$i.txt

      code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o $ofile -d $model --cn_max $cn_max --tree_search $tree_search --init_tree $init_tree --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2  --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5  --verbose $verbose -b 1 > $bdir/std-maxL-"$suffix"-btree-$i
    done
    # Draw the ML tree with bootstapping support
    #Rscript ana/plot-consensus.R $bdir $ofile $dir/results-maxL-tree-sim1-bootstap.pdf
    Rscript ana/plot-trees-all.R -s $bdir -f $mltree -o $dir/results-maxL-tree-"$suffix"-bootstrap.pdf -t "bootstrap" -p "results-maxL-"$suffix"-btree-*txt"
  fi

elif [[ $mode -eq 0 ]]; then
  suffix=sim1-m$model-test-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o $mltree -d $model --cn_max $cn_max --tree_search $tree_search --init_tree $init_tree --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose > $dir/std_svtreeml_"$suffix"
  Rscript ana/plot-trees-all.R -d test/ -b 0 -t "all" -l "xlim" -p "sim-data-*-tree.txt"

elif [[ $mode -eq 2 ]]; then
  # In this mode, the mutation rates have to be specified
  suffix=m$model-"$cons""$maxj"-mode"$mode"-"$prefix"-all"$1"

  code/svtreeml -c $input -t $times -s $Ns --tree_file $tree_file -d $model --cn_max $cn_max --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_svtreeml_"$suffix"

else
  # In this mode, the mutation rates can be specified or not
  suffix=m$model-o"$opt"-"$cons""$maxj"-mode"$mode"-"$prefix"-all"$1"
  mltree=$dir/MaxL-"$prefix"-"$suffix".txt

  code/svtreeml -c $input -t $times -s $Ns --tree_file $tree_file -o $mltree -d $model --cn_max $cn_max --tree_search $tree_search --init_tree $init_tree --dir_itrees $dir_itrees --optim $opt --constrained $cons --fixm $maxj --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode > $dir/std_svtreeml_"$suffix"

  Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim"  #>& /dev/null
fi

# mutree=$dir/results-mu-tree-"$suffix".txt
# code/svtreemu -c $input -t $times -p  $mltree -s $Ns -x $mut -l 0.5  -d $model -o $mutree -m 2000 -r 0.01 > $dir/std_svtreemu_"$suffix"
# #Rscript ana/plot-trees-single.R $dir/results-maxL-mu-tree.txt 0 #>& /dev/null
# Rscript ana/plot-trees-all.R -f $mutree -b 0 -t "single" #>& /dev/null

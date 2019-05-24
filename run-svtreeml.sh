#!/usr/bin/bash

# input parameters
# prefix=sim-data-N"$Ns"-cons"$cons"-model"$model"-1
prefix=sim-data-1
dir=./test
input=$dir/"$prefix"-cn.txt.gz
times=$dir/"$prefix"-rel-times.txt

# number of regions
Ns=5

# Model of evolution. 1: bounded model, 0: JC69 model
model=1
opt=1 # L-BFGS-B method
# opt=0 # simplex method
# Whether or not the tree is constrained by patient age
cons=1
# Whether or not to estimate mutation rate
maxj=1

verbose=0

# seed=111111

# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3

# Specify the mutation rates when they are fixed
mut=0.02
dup_rate=0.00004
del_rate=0.00002

suffix=sim1-m$model-o"$opt"-"$cons""$maxj"

mltree=$dir/results-maxL-tree-"$suffix".txt
code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -d $model -o $mltree -x $mut --dup_rate $dup_rate --del_rate $del_rate --optim $opt --constrained $cons --fixm $maxj --verbose $verbose > $dir/std_svtreeml_"$suffix"
#Rscript ana/plot-trees-single.R $dir/results-maxL-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" -l "xlim"  #>& /dev/null

# mutree=$dir/results-mu-tree-"$suffix".txt
# code/svtreemu -c $input -t $times -p  $mltree -s $Ns -x $mut -l 0.5  -d $model -o $mutree -m 2000 -r 0.01 > $dir/std_svtreemu_"$suffix"
# #Rscript ana/plot-trees-single.R $dir/results-maxL-mu-tree.txt 0 #>& /dev/null
# Rscript ana/plot-trees-all.R -f $mutree -b 0 -t "single" #>& /dev/null


bootstap=1
if [[ $bootstap -eq 1 ]]; then
  # Do bootstapping
  bdir=$dir/bootstrap
  mkdir -p $bdir
  for i in {1..10..1}
  do
    echo $i
    ofile=$bdir/results-maxL-"$suffix"-btree-$i.txt
    code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -d $model -o $ofile -x $mut --dup_rate $dup_rate --del_rate $del_rate --optim $opt --constrained $cons --fixm $maxj --verbose $verbose -b 1 > $bdir/std-maxL-"$suffix"-btree-$i
  done

  # Draw the ML tree with bootstapping support
  #Rscript ana/plot-consensus.R $bdir $ofile $dir/results-maxL-tree-sim1-bootstap.pdf
  Rscript ana/plot-trees-all.R -s $bdir -f $mltree -o $dir/results-maxL-tree-"$suffix"-bootstrap.pdf -t "bootstrap" -p "results-maxL-"$suffix"-btree-*txt"
fi

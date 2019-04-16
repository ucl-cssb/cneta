##!/usr/bin/env bash

# Create a ML tree from copy number profile and do bootstapping

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
minfo=./test/sim-data-1-info.txt
Ns=5
patient=sim1
# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3
model=0 # JC69 model
# opt=1 # L-BFGS-B method
opt=0 # simplex method
cons=0
maxj=0
mut=0.025

bdir=./test/bootstrap
mkdir $bdir
odir=./test

suffix=sim1-m$model-o"$opt"-"$cons""$maxj"

# Create a ML tree
mltree=./test/results-maxL-tree-$suffix.txt
if [ ! -f $mltree ]; then
  code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -d $model -o $mltree -x $mut --optim $opt --constrained $cons --fixm $maxj > ./test/std_svtreeml_"$model"_o"$opt"_"$cons""$maxj"
fi
# > $odir/std-maxL-tree-sim1

# Do bootstapping
for i in {1..10..1}
do
  echo $i
  ofile=$bdir/results-maxL-"$suffix"-btree-$i.txt
  code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o $ofile -b 1 > $bdir/std-maxL-"$suffix"-btree-$i
done

# Draw the ML tree with bootstapping support
#Rscript ana/plot-consensus.R $bdir $ofile $odir/results-maxL-tree-sim1-bootstap.pdf
Rscript ana/plot-trees-all.R -s $bdir -f $ofile -o $odir/results-maxL-tree-"$suffix"-bootstap.pdf -t "bootstrap" -p "results-maxL-"$suffix"-btree-*txt"

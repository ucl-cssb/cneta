##!/usr/bin/env bash

# fit tree model to copy number profile

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

bdir=./test/bootstrap
mkdir $bdir
odir=./test

ofile=$odir/results-maxL-tree-sim1.txt
# code/svtreeml -c $input -t $times -p $Npop -g $Ngen -e $Nstop -o $ofile
#> $odir/std-maxL-tree-sim1

# for i in {1..100..1}
# do
#   echo $i
#   ofile=$bdir/results-maxL-sim1-btree-$i.txt
#   code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -b 1 -o $ofile > $bdir/std-maxL-sim1-btree-$i
# done

Rscript ana/plot-consensus.R $bdir $ofile $odir/results-maxL-tree-sim1-bootstap.pdf

# Find the number of mutations along edges of the simulated tree
less $minfo | sed '1,20d' | grep -v "MUT" | grep -v "eid" | sed '/^$/d' - | awk '{count[$1]++} END{for(e in count) print e, count[e]}' > $odir/sim-data-1-mut-count.txt

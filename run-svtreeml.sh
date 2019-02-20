# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

# evolutionary algorithm
Npop=25
Ngen=50
Nstop=3

code/svtreeml $input $times $Ns $Npop $Ngen $Nstop
Rscript ana/plot-trees.R ./ 0

mv plot-sim-data-1-tree.pdf results-maxL-tree.pdf
mv sim-data-1-tree.txt results-maxL-tree.txt

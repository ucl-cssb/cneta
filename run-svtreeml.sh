# fit tree model to copy number profile

# input parameters
input=./test/sim-data-2-cn.txt.gz
Ns=10

# evolutionary algorithm
Npop=100
Ngen=10000
cnmax=12

code/svtreeml $input $Ns $Npop $Ngen $cnmax
Rscript ana/plot-trees.R ./ 0

mv plot-sim-data-1-tree.pdf results-parsimony-tree.pdf
mv sim-data-1-tree.txt results-parsimony-tree.txt

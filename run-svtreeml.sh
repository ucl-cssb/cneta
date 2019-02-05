# fit tree model to copy number profile

# input parameters
input=./test/sim-data-2-cn.txt.gz
Ns=5

# evolutionary algorithm
Npop=50
Ngen=100
Nstop=5

code/svtreeml $input $Ns $Npop $Ngen $Nstop
Rscript ana/plot-trees.R ./ 0

mv plot-sim-data-1-tree.pdf results-maxL-tree.pdf
mv sim-data-1-tree.txt results-maL-tree.txt

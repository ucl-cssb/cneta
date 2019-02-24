# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3

code/svtreeml $input $times $Ns $Npop $Ngen $Nstop
Rscript ana/plot-trees-single.R results-maxL-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-single.R results-maxL-mu-tree.txt 0 #>& /dev/null



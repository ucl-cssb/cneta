# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3

code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop
Rscript ana/plot-trees-single.R results-maxL-tree.txt 0 #>& /dev/null

mu_est=1.0
code/svtreemu $input $times results-maxL-tree.txt $Ns $mu_est
Rscript ana/plot-trees-single.R results-maxL-mu-tree.txt 0 #>& /dev/null

# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3

#code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o results-maxL-tree.txt
#Rscript ana/plot-trees-single.R results-maxL-tree.txt 0 #>& /dev/null

code/svtreemu -c $input -t $times -p results-maxL-tree.txt -s $Ns -e 1.0 -l 0.5 -o results-maxL-mu-tree.txt
Rscript ana/plot-trees-single.R results-maxL-mu-tree.txt 0 #>& /dev/null

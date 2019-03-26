# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

# evolutionary algorithm
Npop=100
Ngen=50
Nstop=3

code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -o ./test/results-maxL-tree.txt -x 0.025
#Rscript ana/plot-trees-single.R ./test/results-maxL-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-all.R -f ./test/results-maxL-tree.txt -b 0 -t "single" #>& /dev/null

code/svtreemu -c $input -t $times -p ./test/results-maxL-tree.txt -s $Ns -x 0.025 -l 0.5 -o ./test/results-maxL-mu-tree.txt -m 2000 -r 0.01
#Rscript ana/plot-trees-single.R ./test/results-maxL-mu-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-all.R -f ./test/results-maxL-mu-tree.txt -b 0 -t "single" #>& /dev/null

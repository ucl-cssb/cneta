# fit tree model to copy number profile

# input parameters
input=./test/sim-data-1-cn.txt.gz
times=./test/sim-data-1-rel-times.txt
Ns=5

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

suffix=sim1-m$model-o"$opt"-"$cons""$maxj"

mltree=./test/results-maxL-tree-"$suffix".txt
code/svtreeml -c $input -t $times -s $Ns -p $Npop -g $Ngen -e $Nstop -d $model -o $mltree -x $mut --optim $opt --constrained $cons --fixm $maxj > ./test/std_svtreeml_"$suffix"
#Rscript ana/plot-trees-single.R ./test/results-maxL-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-all.R -f $mltree -b 0 -t "single" #>& /dev/null

mutree=./test/results-mu-tree-"$suffix".txt
code/svtreemu -c $input -t $times -p  $mltree -s $Ns -x $mut -l 0.5  -d $model -o $mutree -m 2000 -r 0.01 > ./test/std_svtreemu_"$suffix"
#Rscript ana/plot-trees-single.R ./test/results-maxL-mu-tree.txt 0 #>& /dev/null
Rscript ana/plot-trees-all.R -f $mutree -b 0 -t "single" #>& /dev/null

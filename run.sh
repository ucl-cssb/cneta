#rm sveta; make

#dir="../sim-data/"
dir="./test/"
Nr=4
Ns=10

mkdir $dir
code/sveta $dir $Nr $Ns
Rscript ana/plot-trees.R $dir
Rscript ana/plot-cns.R $dir

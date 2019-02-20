
############################################ tree simulation
# output directory
dir="./test/"

# number of regions. Output will be Nr+1 including germline
# effective population size (increase this for small number of regions)
Ns=5
Ne=2.0
dt=10

# number of samples. Number of tumours to simulate
Nsim=1

# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.5
r2=0.5
r3=0 #0.1
r4=0 #0.1
r5=0 #0.001

# average size of duplications and deletions in bins: L ~ Exp(s)
s1=30
s2=30

mkdir $dir
code/sveta $dir $Ns $Nsim $r1 $r2 $r3 $r4 $r5 $s1 $s2 $Ne $dt
Rscript ana/plot-trees.R $dir 0
Rscript ana/plot-cns.R $dir ana/bin_locations_4401.Rdata


############################################ MLE tree estimate
# input parameters
input="${dir}/sim-data-1-cn.txt.gz"
times="${dir}/sim-data-1-rel-times.txt"

# evolutionary algorithm
Npop=50
Ngen=50
Nstop=5

code/svtreeml $input $times $Ns $Npop $Ngen $Nstop
Rscript ana/plot-trees.R ./ 0

mv plot-sim-data-1-tree.pdf results-maxL-tree.pdf
mv sim-data-1-tree.txt results-maxL-tree.txt


# output directory
dir="./test/"

# number of regions. Output will be Nr+1 including germline
Nr=10

# number of samples. Number of tumours to simulate
Ns=10

# rates of duplication, deletion, chromosome gain, chromosome loss, wgd
r1=0.1
r2=0.1
r3=0.1
r4=0.1
r5=0.001

# average size of duplications and deletions in bins: L ~ Exp(s)
s1=30
s2=30

mkdir $dir
code/sveta $dir $Nr $Ns $r1 $r2 $r3 $r4 $r5 $s1 $s2
Rscript ana/plot-trees.R $dir
Rscript ana/plot-cns.R $dir

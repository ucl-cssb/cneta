
# output directory
dir="./test/"

# number of regions. Output will be Nr+1 including germline
Ns=10

# number of samples. Number of tumours to simulate
Nsim=3

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
code/sveta $dir $Ns $Nsim $r1 $r2 $r3 $r4 $r5 $s1 $s2
Rscript ana/plot-trees.R $dir 1
Rscript ana/plot-cns.R $dir

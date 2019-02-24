
# output directory
dir="./test/"

# number of regions. Output will be Nr+1 including germline
Ns=5

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

# effective population size
Ne=2.0

# observed time step
dt=10

mkdir $dir
code/sveta $dir $Ns $Nsim $r1 $r2 $r3 $r4 $r5 $s1 $s2 $Ne $dt
Rscript ana/plot-trees.R $dir 0 >& /dev/null
Rscript ana/plot-cns.R $dir ana/bin_locations_4401.Rdata >& /dev/null

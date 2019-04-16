

# output directory
dir="./test/"

# number of regions. Output will be Nr+1 including germline
Ns=5

# number of samples. Number of tumours to simulate
Nsim=2

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
code/sveta -o $dir -r $Ns -n $Nsim --dup_rate $r1 --del_rate $r2 --chr_gain $r3 --chr_loss $r4 --wgd $r5 --dup_size $s1 --del_size $s2 -e $Ne -t $dt > ./test/std_test_sveta
#Rscript ana/plot-trees.R $dir 0 >& /dev/null
Rscript ana/plot-trees-all.R -d $dir -b 0 -t "all" # >& /dev/null
# Plot simulated tree with the number of mutations on the branch
Rscript ana/plot-trees-all.R -d $dir -b 1 -t "all" # >& /dev/null
Rscript ana/plot-cns.R -d $dir -b ana/bin_locations_4401.Rdata # >& /dev/null

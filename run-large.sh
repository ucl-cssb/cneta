#!/bin/bash

# take advantage of multiple cores

topdir="./data/"
Nr=4
Ns=10
Nr=10
r1=0.1
r2=0.1
r3=0.1
r4=0.1
r5=0.001
s1=30
s2=30

for i in `seq 1 4`; do
    dir="${topdir}/d${i}"
    mkdir -p $dir 
    code/sveta $dir $Nr $Ns $r1 $r2 $r3 $r4 $r5 $s1 $s2 &
done

wait

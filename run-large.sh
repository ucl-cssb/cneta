#!/bin/bash

topdir="/Users/cbarnes-ihome/data/"
Nr=4
Ns=250

for i in `seq 1 4`; do
    dir="${topdir}/d${i}"
    mkdir $dir 
    ./sveta $dir $Nr $Ns &
done

wait

#!/bin/bash

#PBS -l walltime=1000000000

cd /exports/pierre/EFTofBOSS

python2 MCMC_monopoleonly.py $aa $b $c $d $e $f $g $h $ii

wait


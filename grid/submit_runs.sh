#!/bin/bash
for run in $(seq 0 19)
#for run in 10
do
sbatch --job-name=run$run --output=outfiles/run$run.out mpi_gridcompute.sh $run 20
done


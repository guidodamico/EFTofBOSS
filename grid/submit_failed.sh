#!/bin/bash
for run in $(seq 0 339)
#for run in 10
do
sbatch --job-name=run$run --output=outfiles/failrun$run.out mpi_failed.sh
done


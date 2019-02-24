#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=hns,normal,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --exclude=sh-107-17,sh-101-61,sh-101-62,sh-101-67,sh-105-72, sh-113-15,sh-102-05,sh-103-38,sh-102-01,sh-102-26,sh-107-07
ml load python/2.7.13
ml load py-mpi4py/3.0.0_py27 
ml load gsl/2.3
ml load fftw/3.3.6

mpiexec -n 25 python GridChallenge_runs.py $1 $2


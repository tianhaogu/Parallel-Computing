#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=omp_hw3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --account=eecs587f22_class
#SBATCH --partition=standard

export OMP_NUM_THREADS=4

./hw3_omp.out > submit_4.txt

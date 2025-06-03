#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o /nethome/sdyp11/outputMPI.txt
#SBATCH -e /nethome/sdyp11/errorsMPI.txt
mpirun --bind-to none ejer2 $1 $2 $3 $4
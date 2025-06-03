#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o /nethome/sdyp11/outputS.txt
#SBATCH -e /nethome/sdyp11/errorsS.txt
./sec $1 $2 $3
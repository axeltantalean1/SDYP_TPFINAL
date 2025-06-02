#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o /nethome/sdyp11/outputP.txt
#SBATCH -e /nethome/sdyp11/errorsP.txt
./pth $1 $2 $3 $4
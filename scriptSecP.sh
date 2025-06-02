#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o /nethome/sdyp11/output.txt
#SBATCH -e /nethome/sdyp11/errors.txt
./secPt $1 $2 $3 $4
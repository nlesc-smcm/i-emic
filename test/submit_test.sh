#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 4
#SBATCH --tasks-per-node 16

cd /home/emulder/Projects/I-EMIC/test

srun ./test
#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH -p normal

# link to newest slurm output 
rm -f slurm-latest 
ln -s `ls slurm-* -rt | tail -n 1` slurm-latest  

procs=$1
executable=$2

echo running $executable
echo "#Procs:" $procs

srun -n $procs $executable > dump 

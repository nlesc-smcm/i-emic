#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -n 24

cd ${HOME}/Projects/I-EMIC/rundir_profile

procs=1
echo "Profile on:" ${PLAT} " in:" ${PWD} > profile
# 1,2,4,8,16 procs
for i in {1..5}
do
        echo "Run:   " $i
	echo "#Procs:" $procs
	srun -n $procs ./test > dump
	cat profile_$procs.txt >> profile
	procs=$(($procs*2))
done
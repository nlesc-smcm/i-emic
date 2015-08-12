#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -n 4

cd ${HOME}/Projects/I-EMIC/test

procs=1
echo "Profile on:" ${PLAT} " in:" ${PWD} > profile
# 1,2,4,8,16,32,64
for i in {1..7}
do
        echo $i
	echo $procs
	srun -n $procs ./test > dump
	cat profile_$procs.txt >> profile
	procs=$(($procs*2))
done
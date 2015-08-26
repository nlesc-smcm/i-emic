#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=64

cd ${HOME}/Projects/I-EMIC/rundir_profile

# specify horizontal grid size
n=180
m=80
fname=profile_${n}x${m}
echo writing to $fname
echo $n x $m >> $fname
sed -i "s/Global Grid-Size n.*value.*/Global Grid-Size n\" type=\"int\" value=\"$n\"\/>/" \
	ocean_params.xml
sed -i "s/Global Grid-Size m.*value.*/Global Grid-Size m\" type=\"int\" value=\"$m\"\/>/" \
	ocean_params.xml

procs=1
echo "Profile on:" ${PLAT} " in:" ${PWD} > $fname

# 1,2,4,8,16,32,64 procs
for i in {3..7}
do
	echo "----------------------------------------------------------\
----------------------------------" >> $fname
    echo "Run:   " $i >> $fname
	echo "Run:   " $i 
 	echo "#Procs:" $procs >> $fname
	echo "#Procs:" $procs 
	rm profile_output
	srun -n $procs ./test > dump
	cat profile_output >> $fname
	procs=$(($procs*2))
done

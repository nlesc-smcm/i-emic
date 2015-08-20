#!/bin/bash
#SBATCH -t 60:00:00
#SBATCH -N 1

# Differences between vliegtuig and cartesius ->
#    - srun (cartesius) / mpirun (vliegtuig)
#    - ...

cd ${HOME}/Projects/I-EMIC/rundir_hash

# specify horizontal grid size
n=96
m=96
fname=profile_${n}x${m}
echo writing to $fname
echo $n x $m >> $fname
sed -i "s/Global Grid-Size n.*value.*/Global Grid-Size n\" type=\"int\" value=\"$n\"\/>/" \
	ocean_params.xml
sed -i "s/Global Grid-Size m.*value.*/Global Grid-Size m\" type=\"int\" value=\"$m\"\/>/" \
	ocean_params.xml

procs=8
echo "Test hash on:" ${PLAT} " in:" ${PWD} > $fname

echo "----------------------------------------------------------\
----------------------------------" >> $fname

sed -i "s/Use hashing.*value.*/Use hashing\" type=\"bool\" value=\"false\"\/>/" \
	coupledmodel_params.xml

echo "Run:   " $i >> $fname
echo "Run:   " $i 
echo "#Procs:" $procs >> $fname
echo "#Procs:" $procs
echo "hashing false" >> $fname
echo "hashing false" 
rm profile_output
mpirun -n $procs ./test > dump
cat profile_output >> $fname

echo "----------------------------------------------------------\
----------------------------------" >> $fname

sed -i "s/Use hashing.*value.*/Use hashing\" type=\"bool\" value=\"true\"\/>/" \
	coupledmodel_params.xml

echo "Run:   " $i >> $fname
echo "Run:   " $i 
echo "#Procs:" $procs >> $fname
echo "#Procs:" $procs
echo "hashing true" >> $fname
echo "hashing true" 
rm profile_output
mpirun -n $procs ./test > dump
cat profile_output >> $fname



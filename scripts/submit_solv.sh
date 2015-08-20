#!/bin/bash
#SBATCH -t 60:00:00
#SBATCH -N 1

# Differences between vliegtuig and cartesius ->
#    - srun (cartesius) / mpirun (vliegtuig)
#    - ...

cd ${HOME}/Projects/I-EMIC/rundir_solvertype

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
echo "Test solving scheme on:" ${PLAT} " in:" ${PWD} > $fname

echo "----------------------------------------------------------\
----------------------------------" >> $fname

sed -i "s/Solving scheme.*value.*/Solving scheme\" type=\"char\" value=\"D\"\/>/" \
	coupledmodel_params.xml

echo "Run:   " $i >> $fname
echo "Run:   " $i 
echo "#Procs:" $procs >> $fname
echo "#Procs:" $procs
echo "solving scheme D: decoupled" >> $fname
echo "solving scheme D: decoupled" 
rm profile_output
mpirun -n $procs ./test > dump
cat profile_output >> $fname

echo "----------------------------------------------------------\
----------------------------------" >> $fname

sed -i "s/Solving scheme.*value.*/Solving scheme\" type=\"char\" value=\"E\"\/>/" \
	coupledmodel_params.xml

echo "Run:   " $i >> $fname
echo "Run:   " $i 
echo "#Procs:" $procs >> $fname
echo "#Procs:" $procs
echo "solving scheme E: elimination" >> $fname
echo "solving scheme E: elimination" 
rm profile_output
mpirun -n $procs ./test > dump
cat profile_output >> $fname

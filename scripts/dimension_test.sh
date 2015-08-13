#!/bin/bash
cd ${HOME}/Projects/I-EMIC/rundir_dimension

echo Doubling the horizontal dimensions of the problem... > profile
nprocs=1
for k in 8 16 32 64 128
do
	echo $k x $k >> profile
	sed -i "s/Global Grid-Size n.*value.*/Global Grid-Size n\" type=\"int\" value=\"$k\"\/>/" \
		ocean_params.xml
	sed -i "s/Global Grid-Size m.*value.*/Global Grid-Size m\" type=\"int\" value=\"$k\"\/>/" \
		ocean_params.xml
	time mpirun -np $nprocs ./test > dump
	cat profile_$nprocs.txt >> profile
done

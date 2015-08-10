#!/bin/bash
echo Neumann approximation profile > profile
for k in 2 1 0
do
	sed -i "s/Neumann.*value.*/Neumann approximation\" type=\"int\" value=\"$k\"\/>/" \
		coupledmodel_params.xml
	time mpirun -np 4 ./test > dump
	echo k = $k >> profile
	cat profile_4.txt >> profile
done

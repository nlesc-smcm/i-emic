#!/bin/bash
echo Neumann approximation profile > profile
for k in 0 1 2
do
	sed -i "s/Neumann.*value.*/Neumann approximation\" type=\"int\" value=\"$k\"\/>/" \
		coupledmodel_params.xml
	time ./test > dump
	echo k = $k >> profile
	cat profile_1.txt >> profile
done

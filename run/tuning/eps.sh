#!/bin/bash
for i in 1.1 1.2; 
do  
    sed -i 's/switch steepness.*/switch steepness\" type=\"double\" value=\"'$i'\"\/\>/' seaice_params.xml; 
    while [ `squeue -u emulder | wc -l` -ge 10 ]
    do
     	sleep 1;
    done;
    grep switch seaice_params.xml; 
    bash tuning.sh solarcont/3deg12layers/eps$i; 
done 

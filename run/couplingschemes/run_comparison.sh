#!/bin/bash

# obtain timecode
date=`date +%m%d%y-%H%M`

# original directory
origdir=${PWD}

# set executables
run_coupled=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled
time_coupled=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/time_coupled

if ! [[ -s $run_coupled ]]
then
    echo "executable does not exist: " $run_coupled
    echo "abort"
    exit
fi

#---------------------------------------------------------
# different solving schemes, continuations and transients
ctype=(D Q C)
ptype=(D D B)
clabels=("decoupled" "quasicoupled" "coupledprecB")
tlabels=("time_decoupled" "time_quasicoupled" "time_coupledprecB")

for i in 2 # 1 0
do
    # set combined forcing to zero
    sed -i "s/Combined Forcing.*value.*/Combined Forcing\"    type=\"double\" value=\"0.0\"\/>/" coupledmodel_params.xml

    # set preconditioning type
    sed -i "s/Preconditioning.*value.*/Preconditioning\" type=\"char\" value=\"${ptype[$i]}\"\/>/" coupledmodel_params.xml 

    # set solving scheme type
    sed -i "s/Solving scheme.*value.*/Solving scheme\" type=\"char\" value=\"${ctype[$i]}\"\/>/" coupledmodel_params.xml

    $run_coupled > dump

    bash save.sh 0 ${clabels[$i]} "no description"

    # set combined forcing to continuation destination value (0.1)
    sed -i "s/Combined Forcing.*value.*/Combined Forcing\"    type=\"double\" value=\"0.1\"\/>/" coupledmodel_params.xml

    $time_coupled > dump

    bash save.sh 0 ${tlabels[$i]} "no description"

done

#-----------------------------------------------------------------
# fully coupled continuations, different preconditioning schemes

clabels=("coupledprecD" "coupledprecB" "coupledprecF" "coupledprecG")
ptype=(D B F G)


for i in 0 1 2 3
do
    echo ${ptype[$i]}
    
    # set combined forcing to zero
    sed -i "s/Combined Forcing.*value.*/Combined Forcing\"    type=\"double\" value=\"0.0\"\/>/" coupledmodel_params.xml

    # set preconditioning type
    sed -i "s/Preconditioning.*value.*/Preconditioning\" type=\"char\" value=\"${ptype[$i]}\"\/>/" coupledmodel_params.xml 

    # set solving scheme type (C)
    sed -i "s/Solving scheme.*value.*/Solving scheme\" type=\"char\" value=\"C\"\/>/" coupledmodel_params.xml

    $run_coupled > dump
    bash save.sh 0 ${clabels[$i]} "no description"
done

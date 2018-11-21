#!/bin/bash

# obtain timecode
date=`date +%m%d%y-%H%M`

# original directory
origdir=${PWD}

# set executable
executable=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled

# set number of procs
procs=2

if ! [[ -s $executable ]]
then
    echo "executable does not exist: " $executable
    echo "abort"
    exit
fi

for mdim in 8 16 32 64
do
    # set number of latitude boxes

    for model in ocean atmos seaice coupledmodel
    do
        sed -i "s/Global Grid-Size m.*value.*/Global Grid-Size m\" type=\"int\" value=\"$mdim\"\/>/" ${model}_params.xml
    done

    if [ $# -ge 1 ]
    then
        # create destination dir
        mkdir -pv ${PWD}/$1$mdim
        
        # copy input files to destination dir
        cp -v *.xml $1$mdim
        cp -v *_input.h5 $1$mdim

        # move to destination dir
        cd ${PWD}/$1$mdim
    else
        echo "supply directory"
        echo "abort"
        exit    
    fi

    bash $origdir/setparameters.sh false "Combined Forcing" 1.0 1e-2

    # Run spinup
    if [ -x "$(command -v sbatch)" ]
    then
        bash create_submit.sh "01:00:00" short 1
        sbatch submit.sh $procs $executable
    else        
        mpirun -np $procs $executable > dump
    fi

    for model in ocean atmos seaice
    do
        cp -v ${model}_output.h5 ${model}_input.h5
    done

    bash $origdir/setparameters.sh true "Solar Forcing" 0.0 -1e-2

    # Run solar forcing continuation
    if [ -x "$(command -v sbatch)" ]
    then
        bash create_submit.sh "01:00:00" short 1
        sbatch submit.sh $procs $executable
    else        
        mpirun -np $procs $executable > dump
    fi

    cd $origdir
done

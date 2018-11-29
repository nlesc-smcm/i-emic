#!/bin/bash

# obtain timecode
date=`date +%m%d%y-%H%M`

# original directory
origdir=${PWD}

# set executable
executable=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled

# set number of procs
procs=2

function run
{
    if [ -x "$(command -v sbatch)" ]
    then
	echo "detected slurm, submitting short job"
        bash create_submit.sh short 01:00:00 1
        sbatch submit.sh $1 $2 > jobid
	sed -i "s/.*job //" jobid
	jobid=$(cat jobid)
	sleep 5
	while [ $(squeue -u ${LOGNAME} | grep $jobid | wc -l) -eq 1 ] 
	do
	    echo "job " $jobid " in progress"
	    sleep 10
	done
    else        
	echo "slurm not detected"
        mpirun -np $1 $2 > dump
    fi
}

if ! [[ -s $executable ]]
then
    echo "executable does not exist: " $executable
    echo "abort"
    exit
fi

mdim=(8 16 32)
proc=(2 4 8)

for i in  0 1 2 3
do
    # set number of latitude boxes

    for model in ocean atmosphere seaice coupledmodel
    do
        sed -i "s/Global Grid-Size m.*value.*/Global Grid-Size m\" type=\"int\" value=\"${mdim[$i]}\"\/>/" ${model}_params.xml
    done

    if [ $# -ge 1 ]
    then
        # create destination dir
        mkdir -pv ${PWD}/$1${mdim[$i]}
        
        # copy input files to destination dir
        cp -v *.xml $1${mdim[$i]}
        cp -v *_input.h5 $1${mdim[$i]}

        # move to destination dir
        cd ${PWD}/$1${mdim[$i]}
    else
        echo "supply directory"
        echo "abort"
        exit    
    fi

    bash setparameters.sh false "Combined Forcing" 1.0 1e-2

    # Run spinup
    run ${proc[$i]} $executable

    for model in ocean atmos seaice
    do
        cp -v ${model}_output.h5 ${model}_input.h5
    done

    bash setparameters.sh true "Solar Forcing" 0.0 -1e-2

    # Run solar forcing continuation
    run ${proc[$i]} $executable
    
    cd $origdir
done


#!/bin/bash

# ----------------------------------------------------------------
# Continue from or start run based on xml files in directory <origin> 
#  Put new run in <origin>/<dest>
# ----------------------------------------------------------------
echo "--------------------------------------------"    
if [ $# -lt 1 ]
then
    echo "usage: ./continue <origin>"
    echo "       ./continue <origin> <dest>"
    echo "       ./continue <origin> <dest> <nodetype> <time>"
    echo "       ./continue <origin> <dest> <nodetype> <time> <smartds>"
    echo "        "
    exit
fi

numprocs=8
executable=run_coupled

# figure out some paths
fullexec=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/$executable
scriptsdir=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/scripts\//'`

pwd=${PWD}

if [ $# -ge 2 ]
then
    rdir=$2
else
    rdir=continue
fi

if [ $# -ge 3 ]
then
    type=$3
else
    type=short
fi

if [ $# -ge 4 ]
then
    time=$4
else
    time=01:00:00
fi

if ! [[ -s $1 ]]
then
    echo "No such directory: "$pwd/$1
    exit
else
    echo "Next run can be found in" $1/$rdir
fi

echo $1"/"$rdir >> continue.log

cd $1 
mkdir -p $rdir
mkdir -p $rdir/xml

if ! [[ -s $rdir ]]
then
    echo "Directory creation failed"
    exit
fi    

newds=0
if [[ -s cdata.txt ]] && [ $# -ge 5 ]
then
    # compute new continuation step size
    par0=`tail -n 2 cdata.txt | awk '{print $1}' | head -n 1`
    par1=`tail -n 1 cdata.txt | awk '{print $1}'`
    newds=`awk -v a=$par1 -v b=$par0 'BEGIN{print ((a - b)*100) }'`
    echo "new continuation ds: " $newds
fi

# Put hdf5 files in <dest> directory
if ! [[ -s 'ocean_output.h5' ]]
then
    echo "No hdf5 files in <origin>="$1" to copy!"
    
else
    echo "Copying output hdf5 files in <origin>="$1" to <dest>="$rdir"."
    echo "If loading these files gives problems, try the backup files *_output.h5.bak."
    for i in ocean atmos seaice; 
    do 
        cp  $i'_output.h5' $rdir/$i'_input.h5';
        #cp -v $i'_output.h5.bak' $rdir/$i'_input.h5'; 
    done
fi

# copy xmls twice and cp to continue dir
cp *.xml $rdir
cp *.xml $rdir/xml

cd $rdir

if [ $# -ge 5 ] && [ $newds -ne 0 ] # smart ds
then
    # adjust continuation step size and direction automatically
    echo "Adjust continuation step size and direction automatically based on cdata.txt"
    sed -i "s/initial step size.*value.*/initial step size\"  type=\"double\"  value=\"$newds\"\/>/" continuation_params.xml
fi

if ! [[ -s 'ocean_input.h5' ]]
then
    echo "No input hdf5 files in <dest>="$rdir" directory so we are not loading an existing state!"
    sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"false\"\/>/" ocean_params.xml
    sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"false\"\/>/" coupledmodel_params.xml
else
    echo "Found input hdf5 files in <dest>="$rdir" directory so we are loading an existing state!"        
    sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"true\"\/>/" ocean_params.xml
    sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"true\"\/>/" coupledmodel_params.xml
fi

function run
{    
    if [ -x "$(command -v sbatch)" ]
    then
	    echo "detected slurm, submitting job"
	    echo " type: " $3 " time: " $4 

        bash   $scriptsdir/create_submit.sh $3 $4 1
        sbatch submit.sh $1 $2 > jobid
	    sed -i "s/.*job //" jobid
	    jobid=$(cat jobid)
	    sleep 5
	    if [ $(squeue -u ${LOGNAME} | grep $jobid | wc -l) -eq 1 ] 
	    then
	        echo "job " $jobid " in progress"
	    else
	        echo "job submission failed"
	    fi
    else        
	    echo "slurm not detected"
        mpirun -np $1 $2 > dump
    fi
}

# submit problem (see run() )
run $numprocs $fullexec $type $time
echo "--------------------------------------------"

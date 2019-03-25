#!/bin/bash

# -----------------------------------------
# Continue from run in directory <dir> 
#  Put new run in <dir>/<restart>
# -----------------------------------------

numprocs=16

if [ $# -lt 1 ]
then
    echo "usage: ./continue <dir>"
    echo "       ./continue <dir> <restart>"
    echo "       ./continue <dir> <restart> <type> <time>"
    echo "        "
    exit
fi

function run
{
   
    if [ -x "$(command -v sbatch)" ]
    then
	echo "detected slurm, submitting job"
	echo " type: " $3 " time: " $4 

        bash create_submit.sh $3 $4 1
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

pwd=${PWD}

if [ $# -ge 2 ]
then
    rdir=$2
else
    rdir=restart
fi

if [ $# -ge 3 ]
then
    type=$3
else
    type=normal
fi

if [ $# -ge 4 ]
then
    time=$4
else
    time=5-00:00:00
fi

echo "restart dir is " $dir/$rdir

executable=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled

echo $1"/"$rdir >> continue.log

cd $1 
mkdir -p $rdir

if ! [[ -s $rdir ]]
then
    echo "directory creation failed"
    exit
fi    

# compute new continuation step size
par0=`tail -n 2 cdata.txt | awk '{print $1}' | head -n 1`
par1=`tail -n 1 cdata.txt | awk '{print $1}'`
newds=`awk -v a=$par1 -v b=$par0 'BEGIN{print ((a - b)*100) }'`
echo "new continuation ds: " $newds

# put backup hdf5 files in restart directory
for i in ocean atmos seaice; 
do 
    cp -v $i'_output.h5.bak' $rdir/$i'_input.h5'; 
done

# copy xmls and cp to restart dir
cp -v *.xml $rdir

cd $rdir

# adjust continuation step size and direction
sed -i "s/initial step size.*value.*/initial step size\"  type=\"double\"  value=\"$newds\"\/>/" continuation_params.xml

# just to be sure
sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"true\"\/>/" coupledmodel_params.xml

# submit problem (see run() )
run $numprocs $executable $type $time

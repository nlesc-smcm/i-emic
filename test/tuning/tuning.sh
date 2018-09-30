#!/bin/bash

# #SBATCH --time=96:00:00
#SBATCH --time=1:00:00
#SBATCH --ntasks=24
#SBATCH --nodes=1

# #SBATCH -p fat
# #SBATCH -p normal
#SBATCH -p short

# origin dir
origdir=${PWD}

# This assumes the root dir is called i-emic, otherwise just point to
# the full path here, e.g., <rootdir>/build/src/main/time_coupled
# executable=`echo $origdir | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/time_coupled
executable=`echo $origdir | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled

if ! [[ -s $executable ]]
then
    echo "executable does not exist: " $executable
    echo "abort"
    exit
fi

# cd to original dir
cd ${PWD}

if [ $# -ge 1 ]
then
    # create destination dir
    mkdir -pv ${PWD}/$1
    
    # copy input files to destination dir
    cp -v *.xml $1
    cp -v *_input.h5 $1

    # move to destination dir
    cd ${PWD}/$1
fi

# create filenames
date=`date +%m%d%y-%H%M`
fname=summary_$date
fnameprev=summary_prev_$date
infofile=info_$date
cdatafile=cdata_$date
echo running $executable
echo writing to $logdir/$fname
procs=4

logdir=sbatch_log
mkdir -p $logdir

echo "#Procs:" $procs >> $logdir/$fname
echo "#Procs:" $procs

# link to newest slurm output
rm -f slurm-latest $logdir/slurm-latest
ln -s `ls "$origdir"/slurm-* -rt | tail -n 1` slurm-latest
ln -s `ls "$origdir"/slurm-* -rt | tail -n 1` $logdir/slurm-latest

# backup old data
cat profile_output > $logdir/$fnameprev
cat *.plot >> $logdir/$fnameprev
cp -v info_0.txt $logdir/$infofile
cp -v cdata.txt $logdir/$cdatafile

# backup params
cp -v ocean_preconditioner_params.xml  $logdir/prec_params_$date
cp -v ocean_params.xml                 $logdir/ocean_params_$date
cp -v coupledmodel_params.xml          $logdir/coupledmodel_params_$date
cp -v solver_params.xml                $logdir/solver_params_$date
cp -v atmosphere_params.xml            $logdir/atmosphere_params_$date
cp -v jdqz_params.xml                  $logdir/jdqz_params_$date
cp -v continuation_params.xml          $logdir/continuation_params_$date

# run, second dummy argument triggers different call
if [ $# -ne 2 ]
then
    srun -n $procs $executable > dump
else
    mpirun -np $procs $executable > dump
fi

# when finished the profile is appended to the summary
cat profile_output >> $logdir/$fname

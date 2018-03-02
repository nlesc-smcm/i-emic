#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH -p short

# origin dir
origdir=${PWD}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SHARED_DIR}/i-emic/lib
export PATH=${PATH}:${SHARED_DIR}/i-emic/bin

# create destination dir
mkdir -pv ${HOME}/Projects/i-emic/run/ocean/global/$1

# cd to original dir
cd ${HOME}/Projects/i-emic/run/ocean/global/

# copy input files to destination dir
cp *.xml $1
cp ocean_input.h5 $1

# move to destination dir
cd ${HOME}/Projects/i-emic/run/ocean/global/$1

# create filenames
date=`date +%m%d%y-%H%M`
fname=summary_$date
fnameprev=summary_prev_$date
infofile=info_$date
cdatafile=cdata_$date
echo writing to $logdir/$fname
procs=16

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

# run
srun -n $procs ${HOME}/Projects/i-emic/build/src/main/run_ocean > dump

# when finished the profile is appended to the summary
cat profile_output >> $logdir/$fname

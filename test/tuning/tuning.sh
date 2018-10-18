#!/bin/bash

procs=2          # specify number of processors to use
time=24:00:00    # time (HH:MM:SS) for slurm batch job
#time=01:00:00   
nodes=1          # number of nodes
type=normal      # node type
#type=short      

# append this submit to submit log
echo "tuning.sh " $@ >> tuning.log

# Specify executable.
# This assumes the root dir is called i-emic, otherwise just point to
# the full path here, e.g., <rootdir>/build/src/main/time_coupled

# executable=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/time_coupled
executable=`echo ${PWD} | sed 's/i-emic\/.*/i-emic\/build\/src\//'`main/run_coupled

if ! [[ -s $executable ]]
then
    echo "executable does not exist: " $executable
    echo "abort"
    exit
fi

# cd to original dir
origdir=${PWD}
cd $origdir

if [ $# -ge 1 ]
then
    # create destination dir
    mkdir -pv ${PWD}/$1
    
    # copy input files to destination dir
    cp -v *.xml $1
    cp -v *_input.h5 $1
    cp -v submit.sh $1

    # move to destination dir
    cd ${PWD}/$1
fi

# create filenames
date=`date +%m%d%y-%H%M`
fname=summary_$date
fnameprev=summary_prev_$date
infofile=info_$date
cdatafile=cdata_$date

logdir=sbatch_log
mkdir -p $logdir

echo writing to $logdir/$fname
echo running $executable
echo running $executable >> $logdir/$fname

echo "#Procs:" $procs >> $logdir/$fname
echo "#Procs:" $procs

# backup existing data
cat profile_output > $logdir/$fnameprev
cat *.plot >> $logdir/$fnameprev
cp -v info_0.txt $logdir/$infofile
cp -v cdata.txt $logdir/$cdatafile

# run, second dummy argument triggers ordinary mpirun call
if [ $# -ne 2 ]
then

    # generate submit.sh script
    echo "generating "$origdir/$1"/submit.sh script"

    echo "#!/bin/bash" > submit.sh
    echo "#SBATCH --time="$time >> submit.sh
    echo "#SBATCH --nodes="$nodes >> submit.sh
    echo "#SBATCH -p "$type >> submit.sh
    echo "
# link to newest slurm output 
rm -f slurm-latest 
ln -s \`ls slurm-* -rt | tail -n 1\` slurm-latest  

procs=\$1
executable=\$2

echo running \$executable
echo \"#Procs:\" \$procs

srun -n \$procs \$executable > dump " >> submit.sh
    
    # submit generated script
    sbatch submit.sh $procs $executable
else
    mpirun -np $procs $executable > dump
fi

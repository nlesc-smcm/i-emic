#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --ntasks=64

cd ${HOME}/Projects/i-emic/run/present/coupled/

date=`date +%m%d%y-%H%M`
fname=summary_$date
fnameprev=summary_prev_$date
infofile=info_$date
cdatafile=cdata_$date
profilef=profile_$date
echo writing to $fname

procs=1
echo "Profile on:" ${PLAT} " in:" ${PWD} > $fname

# 1,2,4,8,16,32,64 procs
for i in {3..7}
do
    echo "----------------------------------------------------------\
----------------------------------" >> $fname
    echo "Run:   " $i >> $fname
    echo "Run:   " $i 
    echo "#Procs:" $procs >> $fname
    echo "#Procs:" $procs 
    
    # backup old data
    cat profile_output > $fnameprev
    cat *.plot >> $fnameprev
    cp -v info_0.txt $infofile
    cp -v cdata.txt $cdatafile
    cp -v profile_output $profilef

    rm profile_output
    srun -n $procs run_coupled > dump
    cat profile_output >> $fname
    procs=$(($procs*2))
done

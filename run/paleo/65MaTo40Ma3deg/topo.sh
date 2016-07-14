#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=24
#SBATCH --nodes=1

#SBATCH -p fat

cd ${HOME}/Projects/i-emic/run/paleo/65MaTo40Ma3deg/

date=`date +%m%d%y-%k%M`
fname=summary_$date
fnameprev=summary_prev_$date
infofile=info_$date
echo writing to $fname
procs=24

echo "Global configuration on:" ${PLAT} " in:" ${PWD} > $fname
echo "#Procs:" $procs >> $fname
echo "#Procs:" $procs 

# backup old data
cat profile_output > $fnameprev
cat *.plot >> $fnameprev
cp info_0.txt $infofile

# run
srun -n $procs run_topo > dump

cat profile_output >> $fname

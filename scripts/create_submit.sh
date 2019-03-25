#!/bin/bash

type=$1
time=$2
nodes=$3

# generate submit.sh script
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

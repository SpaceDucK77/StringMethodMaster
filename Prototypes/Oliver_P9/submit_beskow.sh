#!/bin/bash
# Include your allocation number
#SBATCH -A 2019-2-34

# The name of the job in the queue
#SBATCH -J markoStringMethod

# Total number of nodes
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=32

# length in hours
#SBATCH -t 00:59:59

# Receive e-mails when your job starts and ends
#SBATCH --mail-user=oliver.fleetwood@scilifelab.se --mail-type=FAIL

# Output file names for stdout and stderr
#SBATCH -e error.log
#SBATCH -o output.log

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/pdc/vol/anaconda/2019.03/py37/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/pdc/vol/anaconda/2019.03/py37/etc/profile.d/conda.sh" ]; then
        . "/pdc/vol/anaconda/2019.03/py37/etc/profile.d/conda.sh"
    else
       export PATH="/pdc/vol/anaconda/2019.03/py37/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

#module swap PrgEnv-cray PrgEnv-gnu
module load gromacs/2020.1
export OMP_NUM_THREADS=32
export GMX_MAXBACKUP=-1
echo "LOADED MODULES"

rm debug.pickle
rm index.ndx
cp original_index.ndx index.ndx
cmd="conda activate /cfs/klemming/nobackup/o/oliverfl/py37"
echo $cmd
$cmd
cmd="srun -n 1 python VIS_col.py"
#cmd="python VIS_col.py"
echo $cmd
$cmd

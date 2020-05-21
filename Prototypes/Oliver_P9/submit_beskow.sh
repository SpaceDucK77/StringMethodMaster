#!/bin/bash
# Include your allocation number
#SBATCH -A 2019-2-34

# The name of the job in the queue
#SBATCH -J markoStringMethod

# Total number of nodes
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=32

# length in hours
#SBATCH -t 00:09:59

# Receive e-mails when your job starts and ends
#SBATCH --mail-user=oliver.fleetwood@scilifelab.se --mail-type=FAIL

# Output file names for stdout and stderr
#SBATCH -e error.log
#SBATCH -o output.log

module swap PrgEnv-cray PrgEnv-gnu
module load gromacs/2020.1
export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=-1
echo "LOADED MODULES"

cmd="conda activate /cfs/klemming/nobackup/o/oliverfl/py37"
echo $cmd
$cmd
cmd="srun -n 1 python VIS_col.py"
echo $cmd
$cmd

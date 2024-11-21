#!/bin/bash
# =================
# script.sh
# =================

#SBATCH --job-name=boid_num_5
#SBATCH --partition=teach_cpu
#SBATCH --account=PHYS033184
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0:10:00
#SBATCH --mem=100M
# Load modules required for runtime e.g.

module load languages/Intel-OneAPI/2024.0.2


cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=28




# Now run your program with the usual command

echo NON CELL LIST

srun --mpi=pmi2 /user/home/oz21652/Comp401/bc4/bin/5.exe


#!/bin/bash
# =================
# script.sh
# =================

#SBATCH --job-name=boid_num_100
#SBATCH --partition=teach_cpu
#SBATCH --account=PHYS033184
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=1:00:00
#SBATCH --mem=200M
# Load modules required for runtime e.g.

module load languages/Intel-OneAPI/2024.0.2

export OMP_NUM_THREADS=28

cd $SLURM_SUBMIT_DIR


# Now run your program with the usual command

echo NON CELL LIST

srun --mpi=pmi2 /user/home/oz21652/Comp401/bc4/bin/100.exe


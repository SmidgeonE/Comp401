#!/bin/bash
# =================
# script.sh
# =================

#SBATCH --job-name=boid_num_20
#SBATCH --partition=teach_cpu
#SBATCH --account=PHYS033184
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:0:30
#SBATCH --mem=100M
# Load modules required for runtime e.g.

module load languages/Intel-OneAPI/2024.0.2

export OMP_NUM_THREADS=$ {SLURM_CPUS_PER_TASK}

cd $SLURM_SUBMIT_DIR


# Now run your program with the usual command

srun --mpi=pmi2 bc4/bin/20.exe


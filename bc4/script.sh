#!/bin/bash
# =================
# script.sh
# =================

#SBATCH --job-name=test_job
#SBATCH --partition=teach_cpu
#SBATCH --account=PHYS033184
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=7
#SBATCH --time=0:0:10
#SBATCH --mem=100M
# Load modules required for runtime e.g.

module load languages/Intel-OneAPI/2024.0.2

export OMP_NUM_THREADS=$ {SLURM_CPUS_PER_TASK}

cd $SLURM_SUBMIT_DIR


# Now run your program with the usual command

srun --mpi=pmi2 bc4/outputMpi.exe


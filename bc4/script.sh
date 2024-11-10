#!/bin/bash
# =================
# script.sh
# =================

#SBATCH --job-name=test_job
#SBATCH --partition=teach_cpu
#SBATCH --account=PHYS033184
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:0:10
#SBATCH --mem=100M
# Load modules required for runtime e.g.

cd $SLURM_SUBMIT_DIR

# Now run your program with the usual command

./output.exe 1


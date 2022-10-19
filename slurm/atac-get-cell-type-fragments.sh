#!/bin/bash
#SBATCH --job-name="fragments"
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

module load python/3.7.1

scripts/atac-get-cell-type-fragments.py atac atac $SLURM_ARRAY_TASK_ID $slots

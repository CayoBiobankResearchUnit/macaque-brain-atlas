#!/bin/bash
#SBATCH --job-name="fragments"
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

module load python/3.7.1

scripts/atac-get-cell-subtype-fragments.py atac atac $1 $slots

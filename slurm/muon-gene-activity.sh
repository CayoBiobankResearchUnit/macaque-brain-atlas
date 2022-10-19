#!/bin/bash
#SBATCH --job-name="muon"
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/muon-gene-activity.py $1 $2 $SLURM_ARRAY_TASK_ID

#!/bin/bash
#SBATCH --job-name="monocle"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load r/4.0.2

scripts/monocle-recluster-full.R $1 $2 $SLURM_ARRAY_TASK_ID
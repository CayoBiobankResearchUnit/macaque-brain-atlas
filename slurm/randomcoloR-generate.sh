#!/bin/bash
#SBATCH --job-name="randomcoloR"
#SBATCH --time=4:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load r/4.0.2

Rscript scripts/randomcoloR-generate.R $SLURM_ARRAY_TASK_ID
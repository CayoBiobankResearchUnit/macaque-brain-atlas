#!/bin/bash
#SBATCH --job-name="r"
#SBATCH --mem=0
#SBATCH --time=0:30:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load r/4.0.2

scripts/r-recluster-marker-peaks-plot.R $1 $2 $SLURM_ARRAY_TASK_ID

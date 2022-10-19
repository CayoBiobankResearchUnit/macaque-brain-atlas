#!/bin/bash
#SBATCH --job-name="r"
#SBATCH --mem=0
#SBATCH --time=0:30:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load r/4.0.2

scripts/r-recluster-marker-genes-plot.R $1 $2

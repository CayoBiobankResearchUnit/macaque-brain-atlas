#!/bin/bash
#SBATCH --job-name="nnls"
#SBATCH --time=1:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load r/4.0.2

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

scripts/dataset-correlate.R $1 $2 $slots
#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/scanpy-single.py $1 $2 $SLURM_ARRAY_TASK_ID

module purge
module load r/4.0.2

scripts/scanpy-monocle-single.R $1 $2 $SLURM_ARRAY_TASK_ID
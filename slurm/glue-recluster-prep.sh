#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/glue-recluster-preprocess.py $1 $2 $SLURM_ARRAY_TASK_ID
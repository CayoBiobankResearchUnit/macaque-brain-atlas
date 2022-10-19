#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/glue-recluster-transfer-labels-prep.py $1 $2 $3 $SLURM_ARRAY_TASK_ID

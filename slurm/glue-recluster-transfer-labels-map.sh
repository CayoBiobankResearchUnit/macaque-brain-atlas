#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/glue-recluster-transfer-labels-map.py $1 $2 $3 $4 $SLURM_ARRAY_TASK_ID

module purge
module load r/4.0.2
scripts/glue-recluster-summarize-predictions.R $1 $2 $3 $4 $SLURM_ARRAY_TASK_ID

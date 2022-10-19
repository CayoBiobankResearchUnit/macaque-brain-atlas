#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

scripts/glue-recluster-transfer-labels-knn.py $1 $2 $3 $SLURM_ARRAY_TASK_ID $4 $slots

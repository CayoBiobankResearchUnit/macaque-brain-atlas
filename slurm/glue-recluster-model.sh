#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1
module load bedtools/2.30.0

scripts/glue-recluster-model.py $1 $2 $3 $SLURM_ARRAY_TASK_ID
scripts/glue-recluster-regulatory-inference.py $1 $2 $3 $SLURM_ARRAY_TASK_ID

module purge
module load r/4.0.2

scripts/glue-recluster-monocle.R $1 $2 $3 $SLURM_ARRAY_TASK_ID
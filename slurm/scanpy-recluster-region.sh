#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

scripts/scanpy-recluster-region.py $1 $2 $3 $4

module purge
module load r/4.0.2

scripts/monocle-recluster-region.R $1 $2 $3 $4

#scripts/r-recluster-summarize.R $1 $2 $SLURM_ARRAY_TASK_ID
#
#module purge
#module load python/3.7.1
#
#scripts/scanpy-recluster-marker-genes-prep.py $1 $2 $SLURM_ARRAY_TASK_ID

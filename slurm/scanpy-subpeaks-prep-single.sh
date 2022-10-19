#!/bin/bash
#SBATCH --job-name="sc_prep"
#SBATCH --mem=0
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load python/3.7.1

scripts/scanpy-subpeaks-concat.sh $SLURM_ARRAY_TASK_ID

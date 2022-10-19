#!/bin/bash
#SBATCH --job-name="macs3"
#SBATCH --mem=0
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load python/3.7.1

scripts/macs3-call-sub-subpeaks.sh $1 $SLURM_ARRAY_TASK_ID
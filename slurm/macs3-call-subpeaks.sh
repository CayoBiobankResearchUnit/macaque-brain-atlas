#!/bin/bash
#SBATCH --job-name="macs3"
#SBATCH --mem=0
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

scripts/macs3-call-subpeaks.sh $SLURM_ARRAY_TASK_ID

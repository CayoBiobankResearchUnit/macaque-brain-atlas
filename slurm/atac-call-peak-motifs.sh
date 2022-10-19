#!/bin/bash
#SBATCH --job-name="peak-motifs"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

scripts/atac-call-peak-motifs.sh $1 $SLURM_ARRAY_TASK_ID 
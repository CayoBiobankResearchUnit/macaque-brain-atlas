#!/bin/bash
#SBATCH --job-name="cell-motifs"
#SBATCH --mem=100G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load python/3.7.1

scripts/atac-generate-cell-motif-matrix.py $1 atac $SLURM_ARRAY_TASK_ID
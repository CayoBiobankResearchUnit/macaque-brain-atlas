#!/bin/bash
#SBATCH --job-name="seurat"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/seurat-recluster-marker-peaks-prep.py $1 $2 $SLURM_ARRAY_TASK_ID

module purge
module load r/4.0.2
scripts/seurat-recluster-marker-peaks-make.R $1 $2 $SLURM_ARRAY_TASK_ID
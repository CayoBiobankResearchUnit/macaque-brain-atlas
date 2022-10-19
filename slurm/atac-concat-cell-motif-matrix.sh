#!/bin/bash
#SBATCH --job-name="cell-motifs"
#SBATCH --mem=100G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load python/3.7.1

scripts/atac-concat-cell-motif-matrix.py $1 atac

module purge
module load r/4.0.2

scripts/atac-convert-cell-motif-matrix.R $1 atac

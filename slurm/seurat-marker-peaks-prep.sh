#!/bin/bash
#SBATCH --job-name="seurat"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

scripts/seurat-marker-peaks-prep.py $1 $2

module purge
module load r/4.0.2
scripts/seurat-marker-peaks-make.R $1 $2

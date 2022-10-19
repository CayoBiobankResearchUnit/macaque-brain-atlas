#!/bin/bash
#SBATCH --job-name="seurat"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err

module load r/4.0.2

prefix=$1
analysis=$2
ncells=$3

Rscript scripts/monocle-to-seurat.R $prefix $analysis $ncells
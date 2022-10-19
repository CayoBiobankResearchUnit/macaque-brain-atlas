#!/bin/bash
#SBATCH --job-name="monocle"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --exclusive
#SBATCH --mem=0

module load r/4.0.2

prefix=$1

Rscript scripts/atac-qc-clean.R $prefix
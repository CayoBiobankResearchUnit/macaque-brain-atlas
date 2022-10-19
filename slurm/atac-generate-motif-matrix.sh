#!/bin/bash
#SBATCH --job-name="peak-motifs"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00

scripts/atac-generate-motif-matrix.sh $1
#!/bin/bash
#SBATCH --job-name="bbmap"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00

module load bbmap/38.12

scripts/bbmap-index.sh $1
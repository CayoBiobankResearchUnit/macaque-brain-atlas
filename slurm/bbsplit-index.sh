#!/bin/bash
#SBATCH --job-name="bbmap"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load bbmap/38.12

scripts/bbsplit-index.sh $1 $2 $3
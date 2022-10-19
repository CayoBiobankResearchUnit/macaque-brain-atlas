#!/bin/bash
#SBATCH --job-name="nnls"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

scripts/dataset-preprocess.py $1
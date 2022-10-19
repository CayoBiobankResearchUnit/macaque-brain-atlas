#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

# scripts/scanpy-concat.py $1 $2

scripts/scanpy-preprocess.py $1 $2

module purge
module load r/4.0.2

scripts/monocle-preprocess.R $1 $2

#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

scripts/scanpy-marker-peaks-prep.py $1 $2

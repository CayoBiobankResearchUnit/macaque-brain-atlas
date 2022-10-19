#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

scripts/glue-metacells.py $1 $2 $3

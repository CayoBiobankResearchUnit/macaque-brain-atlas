#!/bin/bash
#SBATCH --job-name="glue"
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err

module load python/3.7.1

scripts/glue-transfer-labels-prep.py $1 $2 $3


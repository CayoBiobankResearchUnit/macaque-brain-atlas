#!/bin/bash
#SBATCH --job-name="bbmap"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

id=$(cat data/rna_metadata.txt | tail -n+2 | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f 1)

module load python/3.7.1

scripts/py-bbsplit-parse.py $1 $id $2 $3 $4

module purge
module load r/4.0.2

scripts/r-bbsplit-plot.R $1 $id $2 $3 $4

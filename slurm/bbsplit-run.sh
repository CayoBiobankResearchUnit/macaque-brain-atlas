#!/bin/bash
#SBATCH --job-name="bbmap"
#SBATCH --mem=50G
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

module load bbmap/38.12

id=$(cat data/rna_metadata.txt | tail -n+2 | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f 1)

scripts/bbsplit-run.sh $id $1 $2 $3
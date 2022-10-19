#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load samtools/1.15.1
module load htslib/1.15.1

scripts/samtools-extract-region.sh $SLURM_ARRAY_TASK_ID $1

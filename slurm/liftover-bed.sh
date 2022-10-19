#!/bin/bash
#SBATCH --job-name="liftover"
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load bedtools/2.30.0

bed_prefix=$1

scripts/liftover-bed.sh $bed_prefix $SLURM_ARRAY_TASK_ID

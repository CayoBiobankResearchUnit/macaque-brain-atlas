#!/bin/bash
#SBATCH --job-name="seurat"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load r/4.0.2

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

scripts/seurat-marker-peaks-run.R $1 $2 $SLURM_ARRAY_TASK_ID $slots

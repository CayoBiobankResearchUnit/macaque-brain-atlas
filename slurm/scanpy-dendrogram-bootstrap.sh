#!/bin/bash
#SBATCH --job-name="scanpy"
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load python/3.7.1

prefix=$1
analysis=$2
chunk_size=$3

offset=$((($SLURM_ARRAY_TASK_ID - 1) * $chunk_size))

scripts/scanpy-dendrogram-bootstrap.py $prefix $analysis $chunk_size $offset

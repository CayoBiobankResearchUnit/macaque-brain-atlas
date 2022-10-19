#!/bin/bash
#SBATCH --job-name="sc_prep"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00

module load perl/5.26.0
module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

module load python/3.7.1

parallel -j $slots scripts/scanpy-subpeaks-prep.sh {1} ::: $(seq $1 $2)

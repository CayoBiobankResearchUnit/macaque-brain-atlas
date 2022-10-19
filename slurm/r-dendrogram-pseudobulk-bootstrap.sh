#!/bin/bash
#SBATCH --job-name="sc_prep"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err
#SBATCH --time=7-00:00:00

# module load perl/5.26.0
# module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

module load r/4.0.2

prefix=$1

# parallel -j $slots scripts/r-dendrogram-pseudobulk-bootstrap.R {1} {2} ::: $prefix ::: $(seq 0 1000)

scripts/r-dendrogram-pseudobulk-bootstrap.R $prefix $SLURM_ARRAY_TASK_ID
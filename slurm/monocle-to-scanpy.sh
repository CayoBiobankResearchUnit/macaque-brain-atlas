#!/bin/bash
#SBATCH --job-name="monocle"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

if [ ${1} = "rna" ]; then
module load r/4.0.2
Rscript --min-vsize=1040G --min-nsize=40M scripts/monocle-to-mm.R $1 $2 $SLURM_ARRAY_TASK_ID
module purge
fi

module load python/3.7.1

scripts/scanpy-prep.py $1 $2 $SLURM_ARRAY_TASK_ID

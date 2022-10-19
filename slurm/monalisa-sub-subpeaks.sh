#!/bin/bash
#SBATCH --job-name="monalisa"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err


cell_class=$(cut -f 1 data/monalisa_subpeaks_list.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
cell_cluster=$(cut -f 2 data/monalisa_subpeaks_list.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load r/4.1.0

Rscript scripts/monalisa-sub-subpeaks.R $1 $cell_class $cell_cluster
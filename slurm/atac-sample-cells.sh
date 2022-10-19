#!/bin/bash
#SBATCH --job-name="atac-sample"
#SBATCH --mem=16G
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00

module load r/4.0.2
module load htslib/1.13

this=$(grep '^NSM' data/atac_metadata.txt | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)
cells=$1

Rscript scripts/atac-sample-cells.R $this $cells
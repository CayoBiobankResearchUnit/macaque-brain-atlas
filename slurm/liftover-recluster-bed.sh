#!/bin/bash
#SBATCH --job-name="liftover"
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load bedtools/2.30.0

bed_prefix=$1
atac_prefix=$2

cell_class=$(cut -f 1 data/atacsub_markerpeaks_cell_subclusters.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
cell_cluster=$(cut -f 2 data/atacsub_markerpeaks_cell_subclusters.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

scripts/liftover-recluster-bed.sh $bed_prefix $atac_prefix $cell_class $cell_cluster

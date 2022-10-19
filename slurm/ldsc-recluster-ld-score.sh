#!/bin/bash
#SBATCH --job-name="ldsc"
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

module load perl/5.26.0
module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

bed_prefix=$1

cell_class=$(cut -f 1 data/atacsub_markerpeaks_pass_cell_subclusters.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
cell_cluster=$(cut -f 2 data/atacsub_markerpeaks_pass_cell_subclusters.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load python/3.7.1

parallel -j $slots scripts/ldsc-recluster-ld-score.sh {1} {2} {3} {4} ::: $bed_prefix ::: $cell_class ::: $cell_cluster ::: {1..22}

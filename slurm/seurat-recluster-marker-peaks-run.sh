#!/bin/bash
#SBATCH --job-name="seurat"
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

cell_class=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/atacsub_seurat_cell_subclusters.txt | cut -f 1)
subclass_i=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/atacsub_seurat_cell_subclusters.txt | cut -f 2)

do_parallel=$3

module load r/4.0.2

if [ "$do_parallel" == "1" ]; then
slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)
scripts/seurat-recluster-marker-peaks-run.R $1 $2 $cell_class $subclass_i $slots
else
scripts/seurat-recluster-marker-peaks-run.R $1 $2 $cell_class $subclass_i
fi
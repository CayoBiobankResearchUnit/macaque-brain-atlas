#!/bin/bash

module load python/3.7.1

i=${SLURM_ARRAY_TASK_ID}

source /scratch/nsnyderm/ataac/macs3/bin/activate

macs3 callpeak -f BED -g 2.7e9 \
        -t fragment-files-cells/atac-fragments-*-class${i}.txt.gz \
        -n class${i} --call-summits --nomodel --outdir bed/subpeaks \
        --buffer-size 1000000 --keep-dup all

mv bed/subpeaks/class${i}_peaks.bed bed/subpeaks/merged_peaks-class${i}.bed

# for i in {1..3} {6..13} 15; do
# 	echo class${i}
# 	cp ../../nsnyderm/bbi/atac_peaks/class${i}_peaks.bed.unique bed/subpeaks/unique_peaks-class${i}.bed
# done
# class${i}_peaks.bed.unique
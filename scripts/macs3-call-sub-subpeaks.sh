#!/bin/bash

module load python/3.7.1

i=$1
j=$2

source /scratch/nsnyderm/ataac/macs3/bin/activate

macs3 callpeak -f BED -g 2.7e9 \
        -t fragment-files-clusters/atac-fragments-class${i}_cluster${j}.txt.gz \
        -n class${i}_cluster${j} --call-summits --nomodel --outdir bed/subsubpeaks \
        --buffer-size 1000000 --keep-dup all

cut -f 1-3 bed/subsubpeaks/class${i}_cluster${j}_peaks.narrowPeak | uniq > bed/subsubpeaks/merged_peaks-class${i}_cluster${j}.bed

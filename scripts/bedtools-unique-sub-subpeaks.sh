#!/bin/bash

module load bedtools/2.30.0

i=$1
j=$2

foreground=bed/subsubpeaks/merged_peaks-class${i}_cluster${j}.bed
background=$(ls bed/subsubpeaks/merged_peaks-class${i}_cluster*.bed | grep -v _cluster${j}.bed)

cat $background | bedtools intersect -v -a $foreground -b stdin > bed/subsubpeaks/unique_peaks-class${i}_cluster${j}.bed

# bedtools intersect -wa -a fragment-files-cells/atac-fragments-${this_id}-class${i}.txt.gz


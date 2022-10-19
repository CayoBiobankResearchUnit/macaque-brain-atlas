#!/bin/bash

i=$1
region=$2

# region=18:2923185-3094144

chr=$(printf %02d $(echo $region | cut -d : -f 1))
pos=$(echo $region | cut -d : -f 2 | sed 's/-/_/g')

this_run=$(zcat data/atac_run-data.txt.gz | tail -n+2 | sed -n ${i}p | cut -f 1)
this_id=$(zcat data/atac_run-data.txt.gz | tail -n+2 | sed -n ${i}p | cut -f 2)

mkdir -p bed_region
for i in {1..19}; do
echo class${i}
if [ -f fragment-files-cells/atac-fragments-${this_id}-class${i}.txt.gz ]; then
mkdir -p bed_region/class${i}
mkdir -p bed_region/umi
# echo $region | sed 's/[:-]/'$'\t''/g' | bedtools intersect -wa -a fragment-files-cells/atac-fragments-${this_id}-class${i}.txt.gz -b stdin > bed_region/class${i}/${this_id}_class${i}_region_${chr}_${pos}.bed
zcat fragment-files-cells/atac-fragments-${this_id}-class${i}.txt.gz | wc -l >> bed_region/umi/class${i}_umi.txt
fi
done

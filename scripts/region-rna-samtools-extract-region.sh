#!/bin/bash

i=$1
region=$2

# region=18:2923185-3094144

chr=$(printf %02d $(echo $region | cut -d : -f 1))
pos=$(echo $region | cut -d : -f 2 | sed 's/-/_/g')

this_run=$(zcat data/rna_run-data.txt.gz | tail -n+2 | sed -n ${i}p | cut -f 1)
this_id=$(zcat data/rna_run-data.txt.gz | tail -n+2 | sed -n ${i}p | cut -f 2)

if [ $this_run = "Snyder-Mackler_RNA3-019" ]; then
bam_path=/data/CEM/smacklab/data/bbi/Snyder-Mackler_RNA3-019-nova_data-rerun/bams
else
bam_path=/data/CEM/smacklab/data/bbi/${this_run}-nova_data/bams
fi

if [ $(echo ${this_id#*NSM} | sed 's/\.[a-z]$//g') -lt 100 ]; then
bam_file=NSM${this_id#*NSM0}.bam
else
bam_file=${this_id}.bam
fi

# samtools sort -o bam/${this_id}.sort.bam ${bam_path}/${bam_file}

cp ${bam_path}/${bam_file} bam/${this_id}.bam
samtools index -b bam/${this_id}.bam

mkdir -p bam_region
samtools view -b bam/${this_id}.bam ${region} > bam_region/${this_id}_region_${chr}_${pos}.bam

cat \
<(samtools view -H bam_region/${this_id}_region_${chr}_${pos}.bam) \
<(paste \
<(samtools view bam_region/${this_id}_region_${chr}_${pos}.bam | cut -f 1 | cut -d '|' -f 3-5 | sed 's/|/_/g' | sed 's/^/'${this_id}'-/g') \
<(samtools view bam_region/${this_id}_region_${chr}_${pos}.bam | cut -f 2-15)) | bgzip -c > bam_region/${this_id}_renamed_region_${chr}_${pos}.sam

samtools view -o bam_region/${this_id}_renamed_region_${chr}_${pos}.bam -b bam_region/${this_id}_renamed_region_${chr}_${pos}.sam

samtools index bam_region/${this_id}_renamed_region_${chr}_${pos}.bam

# mkdir -p cell_lists
# for i in {1..19}; do
# echo class${i}
# this_cell=$(sed -n ${i}p stats/clusters/rna-final-cellclasses-levels.txt)
# cat stats/clusters/rna-final-cellclasses.txt | grep "${this_cell}$" | cut -f 1 > cell_lists/class${i}.txt
# done

for i in {1..19}; do
echo class${i}
mkdir -p bam_region/class${i}
samtools view -b -o bam_region/class${i}/${this_id}_class${i}_region_${chr}_${pos}.bam -N cell_lists/class${i}.txt bam_region/${this_id}_renamed_region_${chr}_${pos}.bam
done

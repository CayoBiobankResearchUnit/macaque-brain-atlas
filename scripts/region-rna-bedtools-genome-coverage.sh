#!/bin/bash

region=$1

# region=18:2923185-3094144

chr=$(printf %02d $(echo $region | cut -d : -f 1))
pos=$(echo $region | cut -d : -f 2 | sed 's/-/_/g')

for i in {1..19}; do
echo class${i}
this_cell=$(sed -n ${i}p stats/clusters/rna-final-cellclasses-levels.txt)
bam_list=$(ls bam_region/class${i}/*.bam)
mkdir -p bam_cell
samtools merge -f -o bam_cell/class${i}_region_${chr}_${pos}.bam $bam_list
samtools sort -o bam_cell/class${i}_region_${chr}_${pos}.sorted.bam bam_cell/class${i}_region_${chr}_${pos}.bam
done

# Calculate regions
pos1=$(echo $region | cut -d ':' -f 2 | cut -d '-' -f 1)
pos2=$(echo $region | cut -d ':' -f 2 | cut -d '-' -f 2)

bin_size=100

window1=$(((($pos1+$(($bin_size-1))) / $bin_size) * $bin_size))
window2=$((($pos2/ $bin_size) * $bin_size))

mkdir -p regions
rm -rf regions/${chr}_${pos}_bins.bed
for i in $(seq $window1 $bin_size $(($window2-$bin_size))); do
echo ${chr}$'\t'${i}$'\t'$(($i + $bin_size)) >> regions/${chr}_${pos}_bins.bed
done


# Bedtools

mkdir -p coverage
for i in {1..19}; do
# echo class${i}
this_cell=$(sed -n ${i}p stats/clusters/rna-final-cellclasses-levels.txt)
echo $this_cell
bedtools coverage -a regions/${chr}_${pos}_bins.bed -b bam_cell/class${i}_region_${chr}_${pos}.sorted.bam > coverage/class${i}_${chr}_${pos}_bins.cov.txt
done

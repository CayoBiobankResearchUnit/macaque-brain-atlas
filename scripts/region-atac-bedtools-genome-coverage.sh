#!/bin/bash

region=$1

# region=18:2923185-3094144

chr=$(printf %02d $(echo $region | cut -d : -f 1))
pos=$(echo $region | cut -d : -f 2 | sed 's/-/_/g')

mkdir -p bed_cell
for i in {1..19}; do
echo class${i}
this_cell=$(sed -n ${i}p stats/clusters/rna-final-cellclasses-levels.txt)
if [ -f bed_region/class${i}/NSM*_class${i}_region_${chr}_${pos}.bed ]; then
cat bed_region/class${i}/NSM*_class${i}_region_${chr}_${pos}.bed | sort -k 1,1 -k2,2n > bed_cell/class${i}_region_${chr}_${pos}.bed
fi
done

# # Calculate regions
# pos1=$(echo $region | cut -d ':' -f 2 | cut -d '-' -f 1)
# pos2=$(echo $region | cut -d ':' -f 2 | cut -d '-' -f 2)
# 
# bin_size=100
# 
# window1=$(((($pos1+$(($bin_size-1))) / $bin_size) * $bin_size))
# window2=$((($pos2/ $bin_size) * $bin_size))
# 
# mkdir -p regions
# rm -rf regions/${chr}_${pos}_bins.bed
# for i in $(seq $window1 $bin_size $(($window2-$bin_size))); do
# echo ${chr}$'\t'${i}$'\t'$(($i + $bin_size)) >> regions/${chr}_${pos}_bins.bed
# done


# Bedtools

mkdir -p coverage
for i in {1..19}; do
# echo class${i}
this_cell=$(sed -n ${i}p stats/clusters/rna-final-cellclasses-levels.txt)
echo $this_cell
if [ -f bed_cell/class${i}_region_${chr}_${pos}.bed ]; then
bedtools coverage -a regions/${chr}_${pos}_bins.bed -b bed_cell/class${i}_region_${chr}_${pos}.bed > coverage/atac_class${i}_${chr}_${pos}_bins.cov.txt
fi
done

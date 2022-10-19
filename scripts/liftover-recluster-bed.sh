#!/bin/bash
# 🐛

root=$(pwd)
bed_prefix=$1
atac_prefix=$2
cell_class=$3
cell_cluster=$4

# bed_prefix=marker_peaks

if [ $bed_prefix = 'marker_peaks' ]; then
bedFile=stats/markerpeaks/bed/${atac_prefix}_class${cell_class}_cluster${cell_cluster}_marker_peaks.bed
elif [ \( $bed_prefix = 'merged_peaks' \) -o \( $bed_prefix = 'unique_peaks' \) -o \( $bed_prefix = 'merged_unique_peaks' \) -o \( $bed_prefix = 'top_1p_peaks' \) ]; then
bedFile=bed/subsubpeaks/${bed_prefix}-class${cell_class}_cluster${cell_cluster}.bed
fi

if [ $(grep -c '^chr' $bedFile | xargs) -eq 0 ]; then
	mkdir -p bed/liftover
	cat $bedFile | sed 's/^/chr/g' > bed/liftover/${bed_prefix}_class${cell_class}_cluster${cell_cluster}_liftover.bed
	inBed=bed/liftover/${bed_prefix}_class${cell_class}_cluster${cell_cluster}_liftover.bed
else
	inBed=$bedFile
fi

mkdir -p bed/fail
mkdir -p bed/success

bin/liftOver \
	$inBed \
	data/rheMac10ToHg19.over.chain.gz \
	bed/success/${bed_prefix}_hg19_success_class${cell_class}_cluster${cell_cluster}.bed \
	bed/fail/${bed_prefix}_hg19_fail_class${cell_class}_cluster${cell_cluster}.bed \
	-bedPlus=2 -tab

# # Fix patched chromosome names
# sed 's/\(chr[0-9XY]*\)_[a-z]\{2\}[0-9]\{6\}_fix'$'\t''/\1'$'\t''/g' \
# 	bed/success/${bed_prefix}_hg19_success_class${cell_class}.bed | \
# 	bedtools sort > bed/${bed_prefix}_hg19_class${cell_class}.bed

# Drop peaks in patched scaffolds
grep -v 'chr[0-9XY]*_[a-z]\{2\}[0-9]\{6\}_fix'$'\t' \
	bed/success/${bed_prefix}_hg19_success_class${cell_class}_cluster${cell_cluster}.bed | \
	bedtools sort > bed/${bed_prefix}_hg19_class${cell_class}_cluster${cell_cluster}.bed

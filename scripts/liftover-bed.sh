#!/bin/bash
# 🐛

root=$(pwd)
bed_prefix=$1
cell_class=$2

# bed_prefix=marker_peaks

if [ $(grep -c '^chr' bed/${bed_prefix}_class${cell_class}.bed | xargs) -eq 0 ]; then
	mkdir -p bed/liftover
	cat bed/${bed_prefix}_class${cell_class}.bed | sed 's/^/chr/g' > bed/liftover/${bed_prefix}_class${cell_class}_liftover.bed
	inBed=bed/liftover/${bed_prefix}_class${cell_class}_liftover.bed
else
	inBed=bed/${bed_prefix}_class${cell_class}.bed
fi

mkdir -p bed/fail
mkdir -p bed/success

bin/liftOver \
	$inBed \
	data/rheMac10ToHg19.over.chain.gz \
	bed/success/${bed_prefix}_hg19_success_class${cell_class}.bed \
	bed/fail/${bed_prefix}_hg19_fail_class${cell_class}.bed \
	-bedPlus=2 -tab

# # Fix patched chromosome names
# sed 's/\(chr[0-9XY]*\)_[a-z]\{2\}[0-9]\{6\}_fix'$'\t''/\1'$'\t''/g' \
# 	bed/success/${bed_prefix}_hg19_success_class${cell_class}.bed | \
# 	bedtools sort > bed/${bed_prefix}_hg19_class${cell_class}.bed

# Drop peaks in patched scaffolds
grep -v 'chr[0-9XY]*_[a-z]\{2\}[0-9]\{6\}_fix'$'\t' \
	bed/success/${bed_prefix}_hg19_success_class${cell_class}.bed | \
	bedtools sort > bed/${bed_prefix}_hg19_class${cell_class}.bed
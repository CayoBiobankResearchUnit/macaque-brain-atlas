#!/bin/bash
# 🐛

# Documentation for ldsc: https://github.com/bulik/ldsc
source activate ldsc

bed_prefix=$1
cell_class=$2
cell_cluster=$3
chrom=$4

# Make an output folder
mkdir -p stats/ldsc/${bed_prefix}

ldsc/make_annot.py \
	--bed-file bed/${bed_prefix}_hg19_class${cell_class}_cluster${cell_cluster}.bed \
	--bimfile data/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
	--annot-file stats/ldsc/${bed_prefix}/class${cell_class}_cluster${cell_cluster}.${chrom}.annot.gz

ldsc/ldsc.py \
	--l2 \
	--bfile data/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
	--print-snps data/1000G/1000G_snps/1000G.${chrom}.snp \
	--ld-wind-cm 1 \
	--annot stats/ldsc/${bed_prefix}/class${cell_class}_cluster${cell_cluster}.${chrom}.annot.gz \
	--thin-annot \
	--out stats/ldsc/${bed_prefix}/class${cell_class}_cluster${cell_cluster}.${chrom}

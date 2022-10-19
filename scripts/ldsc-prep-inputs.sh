#!/bin/bash
# 🐛

root=$(pwd)

# Follow this as guide
# https://atlas.gs.washington.edu/mouse-atac/docs/#ld-score-regression

# Scripts here:
# https://github.com/shendurelab/mouse-atac/tree/master/ldscore_regression

inMergedPeaks=${root}/data/merged_peaks.bed

# if [ -d ${root}/shendurelab/mouse-atac ]; then
# 	mkdir shendurelab
# 	cd shendurelab
# 	git clone https://github.com/shendurelab/mouse-atac.git
# 	cd ${root}
# fi
# 
# cd ${root}/shendurelab/mouse-atac

# Obsolete, but keep just in case
if [ ! -f ${root}/data/hg19ToRheMac10.over.chain.gz ]; then
	wget -O ${root}/data/hg19ToRheMac10.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToRheMac10.over.chain.gz
fi
# This below is the liftover chain we'll use
if [ ! -f ${root}/data/rheMac10ToHg19.over.chain.gz ]; then
	wget -O ${root}/data/rheMac10ToHg19.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/rheMac10/liftOver/rheMac10ToHg19.over.chain.gz
fi

# Grab liftover
if [ ! -f ${root}/bin/liftOver ]; then
	mkdir -p ${root}/bin
	wget -O ${root}/bin/liftOver https://github.com/ENCODE-DCC/kentUtils/raw/master/bin/linux.x86_64/liftOver
	chmod +x ${root}/bin/liftOver
fi

# 1000G inputs (https://alkesgroup.broadinstitute.org/LDSCORE)

# Make 1000G inputs

if [ ! -d ${root}/data/1000G ]; then
	mkdir ${root}/data/1000G
fi

inputs="\
1000G_Phase3_baseline_v1.2_ldscores.tgz \
1000G_Phase3_baselineLD_v2.2_ldscores.tgz \
1000G_Phase3_cell_type_groups.tgz \
1000G_Phase3_plinkfiles.tgz \
hapmap3_snps.tgz \
1000G_Phase3_weights_hm3_no_MHC.tgz \
1000G_Phase3_frq.tgz \
"

base_url=https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE

for i in $inputs; do
if [ ! -f ${root}/data/1000G/${i} ]; then
	wget -O ${root}/data/1000G/${i} ${base_url}/${i}
	cd ${root}/data/1000G/
	if [ $i == 1000G_Phase3_baselineLD_v2.2_ldscores.tgz ]; then
		mkdir -p ${root}/data/1000G/baselineLD_v2.2
		mv ${root}/data/1000G/${i} ${root}/data/1000G/baselineLD_v2.2/${i}
		cd ${root}/data/1000G/baselineLD_v2.2/
		tar -xvzf ${i}
		mv ${i} ${root}/data/1000G/
	else
		tar -xvzf ${i}
	fi
	cd ${root}
fi
done

# Cell type groups
# file_num	cell_type
# 1  Adrenal_Pancreas.bed
# 2  Cardiovascular.bed
# 3  CNS.bed
# 4  Connective_Bone.bed
# 5  GI.bed
# 6  Hematopoietic.bed
# 7  Kidney.bed
# 8  Liver.bed
# 9  Other.bed
# 10 SkeletalMuscle.bed

if [ ! -f ${root}/data/1000G/sumstats_full.txt ]; then
	if [ -f ${root}/data/1000G/all_sumstats.txt ]; then
		sumstats_base_url=https://alkesgroup.broadinstitute.org/LDSCORE/all_sumstats
		# Edit file data/1000G/all_sumstats.txt with all files in https://alkesgroup.broadinstitute.org/LDSCORE/all_sumstats/
		mkdir -p ${root}/data/1000G/all_sumstats
		for i in $(cat ${root}/data/1000G/all_sumstats.txt); do
		wget -O ${root}/data/1000G/all_sumstats/${i} ${sumstats_base_url}/${i}
		done

		cat <(echo phenotype$'\t'sumstats) \
		<(paste \
			<(ls -1 ${root}/data/1000G/all_sumstats | sed 's/\.sumstats$//g') \
			<(ls -1 ${root}/data/1000G/all_sumstats | sed 's:\(.*\):'${root}'/data/1000G/all_sumstats/\1:g')) > \
		${root}/data/1000G/sumstats_full.txt
		echo 'Manually filter file "'${root}/data/1000G/sumstats_full.txt'" to create sumstats files.'
	else
		echo 'File '${root}/data/1000G/all_sumstats.txt' must exist'
	fi
elif [ -f ${root}/data/1000G/sumstats_liberal.txt ]; then
	# If data/1000G/sumstats_liberal.txt exists, split into chunks
	mkdir -p ${root}/data/1000G/sumstats_individual
	for i in $(seq 1 $(($(wc -l ${root}/data/1000G/sumstats_liberal.txt | cut -d ' ' -f 1 | xargs) - 1))); do
	echo $i
	cat <(head -n1 ${root}/data/1000G/sumstats_liberal.txt) <(sed -n $(($i + 1))p ${root}/data/1000G/sumstats_liberal.txt) > ${root}/data/1000G/sumstats_individual/sumstats-$(sed -n $(($i + 1))p ${root}/data/1000G/sumstats_liberal.txt | cut -f 1).txt
	done
fi

# Extract SNP lists from baseline files
if [ ! -d data/1000G/1000G_snps ]; then
mkdir -p data/1000G/1000G_snps
for i in {1..22}; do
echo $i
zcat data/1000G/baselineLD_v2.2/baselineLD.${i}.l2.ldscore.gz | cut -f 2 | tail -n+2 > data/1000G/1000G_snps/1000G.${i}.snp
done
fi

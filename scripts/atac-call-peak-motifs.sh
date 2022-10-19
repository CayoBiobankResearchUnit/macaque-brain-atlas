#!/bin/bash
# 🐛

root=$(pwd)
# this=$(grep '^NSM' data/atac_metadata.txt | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)
prefix=$1
i=$(($2-1))

module load bedtools2/2.24.0

inMergedPeaks=${root}/data/merged_peaks.bed

cd ${root}/bbi/bbi-sciatac-analyze

pipeline_path=$(pwd)

module unload python
PYTHONPATH=''
module load gcc/8.2.0
module load python/3.7.1
module load htslib/1.13

activate

script_dir=src

fasta=${root}/data/Macaca_mulatta.Mmul_10.dna.toplevel.fa.finished

# if [ ! -f ${root}/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt ]; then
# 	wget -O ${root}/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt
# fi
# if [ ! -f ${root}/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.pfm ]; then
# 	cat ${root}/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt | sed 's/^[ACGT] *\[ *\([0-9 ]*\) \]/\1/g' | sed 's/ \{1,\}/\t/g' > ${root}/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.pfm
# fi

motifs=${root}/data/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm
mkdir -p ${root}/peak_motifs

python ${script_dir}/call_peak_motifs.py \
	${fasta} \
	${inMergedPeaks} \
	${motifs} \
	${root}/peak_motifs/atac-gc_$(printf "%02d\n" $(($i+1)))-peak_calls.bb  \
	--gc_bin $i \
	--pwm_threshold 1e-7


#!/bin/bash
# 🐛

root=$(pwd)
prefix=$1

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

inPeakCalls=$(/usr/bin/ls ${root}/peak_motifs/${prefix}-gc_*-peak_calls.bb | xargs)

outPeakTfMatrix=${root}/sparse-matrices/${prefix}-peak_tf.mtx.gz

python ${script_dir}/generate_motif_matrix.py \
	--peak_motif_files ${inPeakCalls} \
	--fasta ${fasta} \
	--peaks ${inMergedPeaks} \
	--motifs ${motifs} \
	--peak_tf_matrix ${outPeakTfMatrix}



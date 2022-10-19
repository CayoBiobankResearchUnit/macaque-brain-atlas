#!/bin/bash

root=$(pwd)
this=$(grep '^NSM' data/atac_metadata.txt | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

inTranspositionSites=$(/usr/bin/ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data/NSM{201..310}/get_unique_fragments/*-transposition_sites.bed.gz | grep $this | xargs)
inCellWhitelist=$(/usr/bin/ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data/NSM{201..310}/call_cells/*-called_cells_whitelist.txt | grep $this | xargs)

# These files have to be generated
# Merged peaks bed file
# ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0[01][01279]*nova_data/*/call_peaks/*-merged_peaks.bed

if [ ! -f ${root}/data/merged_peaks.bed ]; then
	for i in $(echo /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data | xargs); do
	cat ${i}/*/call_peaks/*-merged_peaks.bed; done | \
		cut -f1-3 \ | sort -k1,1V -k2,2n -k3,3n | bedtools merge -i - | sort -k1,1V -k2,2n -k3,3n >  merged_peaks.bed
fi

inMergedPeaks=${root}/data/merged_peaks.bed

# head -n 22 ${root}/../cayo_bulk_brain/genomes/Mmul_10.dna.fa.fai | cut -f 1-2 > ${root}/data/mmul_sizes.txt
inChromosomeSizes=${root}/data/mmul_sizes.txt

# git clone https://github.com/bbi/bbi-sciatac-analyze.git

cd bbi/bbi-sciatac-analyze

pipeline_path=$(pwd)

# Had to update the modules below
# source ${pipeline_path}/load_python_env_reqs.sh

module unload python
PYTHONPATH=''
module load gcc/8.2.0
module load python/3.7.1
module load htslib/1.13

# source ${script_dir}/python_env/bin/activate
activate

mkdir -p ${root}/sparse-matrices

outPeakMatrix=${root}/sparse-matrices/${this}-peak_matrix.mtx.gz

script_dir=src

# Need to add sample IDs, do one sample for now
# cat /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0[01][01279]*nova_data/*/get_unique_fragments/*-transposition_sites.bed.gz /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0[01][1]*nova_data/*/get_unique_fragments/NSM2*-transposition_sites.bed.gz

module load bedtools2/2.24.0

python ${script_dir}/generate_sparse_matrix.py \
--transposition_sites_intersect <(bedtools intersect -sorted -g ${inChromosomeSizes} -a ${inMergedPeaks} -b ${inTranspositionSites} -wa -wb) \
--intervals ${inMergedPeaks} \
--cell_whitelist ${inCellWhitelist} \
--matrix_output ${outPeakMatrix}

mv ${root}/sparse-matrices/${this}-peak_matrix.columns.txt ${root}/sparse-matrices/.${this}-peak_matrix.columns.txt

awk '{new_var="'$this'-"$1; print new_var}' ${root}/sparse-matrices/.${this}-peak_matrix.columns.txt > ${root}/sparse-matrices/${this}-peak_matrix.columns.txt

# Transfer count reports
inCountReport=$(/usr/bin/ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data/*/count_report/*-count_report.txt | grep $this | xargs)

mkdir -p count-reports

awk '{if (NR==1) printf $0"\n"; else printf $0="'${this}'-"$0"\n"}' $inCountReport > count-reports/${this}-count_report.txt

mkdir -p fragment-files
# Transfer fragment files
inFragmentFile=$(/usr/bin/ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data/*/get_unique_fragments/*-fragments.txt.gz | grep $this | xargs)

zcat $inFragmentFile | awk '{print $1"\t"$2"\t"$3"\t'${this}-'"$4"\t"$5}' | bgzip -c > fragment-files/${this}-fragments.txt.gz
tabix -p bed fragment-files/${this}-fragments.txt.gz

mkdir -p transposition-sites
# Transfer transposition
inTranspositionSites=$(/usr/bin/ls /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0{07,{09..12}}_nova_data/*/get_unique_fragments/*-transposition_sites.bed.gz | grep $this | xargs)

zcat $inTranspositionSites | awk '{print $1"\t"$2"\t"$3"\t'${this}-'"$4}' | bgzip -c > transposition-sites/${this}-transposition_sites.bed.gz
tabix -p bed transposition-sites/${this}-transposition_sites.bed.gz

# Create copies of atac data to match syntax of rna data
zcat $outPeakMatrix > mm/atac-${this}-expr.mm
cat ${root}/sparse-matrices/${this}-peak_matrix.columns.txt | gzip -c > mm/atac-${this}-expr.cols.txt.gz
paste -d '_' <(cut -f 1 data/merged_peaks.bed) <(cut -f 2 data/merged_peaks.bed) <(cut -f 3 data/merged_peaks.bed) | gzip -c > mm/atac-${this}-expr.rows.txt.gz



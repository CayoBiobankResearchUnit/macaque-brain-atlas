#!/bin/bash

root=$(pwd)
this_job=${1}

this=$(cat data/transposition_file_sizes.txt | sed -n ${this_job}p | cut -f 2)
cell_class=$(cat data/transposition_file_sizes.txt | sed -n ${this_job}p | cut -f 3)

inTranspositionSites=${root}/transposition-sites-cells/atac-transposition-sites-${this}-class${cell_class}.bed.gz
inCellWhitelist=${root}/whitelists-cells/atac-whitelists-${this}-class${cell_class}.txt

inMergedPeaks=${root}/bed/subpeaks/merged_peaks-class${cell_class}.bed

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

mkdir -p ${root}/sparse-matrices/class${cell_class}-matrix

outPeakMatrix=${root}/sparse-matrices/class${cell_class}-matrix/${this}-class${cell_class}-peak_matrix.mtx.gz

script_dir=src

# Need to add sample IDs, do one sample for now
# cat /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0[01][01279]*nova_data/*/get_unique_fragments/*-transposition_sites.bed.gz /data/CEM/smacklab/data/bbi/Snyder-Mackler_ATAC3-0[01][1]*nova_data/*/get_unique_fragments/NSM2*-transposition_sites.bed.gz

module load bedtools2/2.24.0

rm -rf ${outPeakMatrix}
python ${script_dir}/generate_sparse_matrix.py \
--transposition_sites_intersect <(bedtools intersect -sorted -g ${inChromosomeSizes} -a ${inMergedPeaks} -b ${inTranspositionSites} -wa -wb) \
--intervals ${inMergedPeaks} \
--cell_whitelist ${inCellWhitelist} \
--matrix_output ${outPeakMatrix}

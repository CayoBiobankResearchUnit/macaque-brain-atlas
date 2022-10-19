#!/bin/bash
# 🐛

# Documentation for ldsc: https://github.com/bulik/ldsc
source activate ldsc

bed_prefix=$1
cell_class=$2
phenotype=$3

# Make an output folder
mkdir -p stats/ldsc/${bed_prefix}/scores

sumstats_file=sumstats_liberal.txt

phenotype_code=$(sed -n $(($phenotype + 1))p data/1000G/${sumstats_file} | cut -f 1 | xargs)
phenotype_file=$(sed -n $(($phenotype + 1))p data/1000G/${sumstats_file} | cut -f 2 | xargs)

baseline_folder=baselineLD_v2.2
baseline_prefix=baselineLD.

ldsc/ldsc.py \
	--h2 ${phenotype_file} \
	--w-ld-chr data/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--ref-ld-chr stats/ldsc/${bed_prefix}/class${cell_class}.,data/1000G/${baseline_folder}/${baseline_prefix} \
	--frqfile-chr data/1000G/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out stats/ldsc/${bed_prefix}/scores/class${cell_class}-${phenotype_code}

# Write summary stats
python scripts/ldsc-partitioned-h2-summarize.py $bed_prefix $cell_class $phenotype

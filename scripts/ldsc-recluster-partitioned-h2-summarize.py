#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import os
import re

arguments = sys.argv
# arguments = ['ldsc-partitioned-h2-summarize.py','1','1']

bed_prefix = arguments[1]
cell_class = int(arguments[2])
cell_cluster = int(arguments[3])
phenotype = int(arguments[4])

sumstats_file = os.path.join('data/1000G','sumstats_liberal.txt')

sumstats = pd.read_csv(sumstats_file,header=0,index_col=None,sep='\t')

# Turn phenotype index into phenotype label
phenotype_code = sumstats.loc[phenotype-1]['phenotype']

# Turn class index into class label
cell_class_name = pd.read_csv(os.path.join('stats/clusters','rna-final-cellclasses-levels.txt'),index_col=None,header=None,sep='\t')[0].to_list()[cell_class-1]

# Turn cluster number into label
cell_cluster_name = cell_class_name + ' ' + str(cell_cluster)

results_file = os.path.join('stats/ldsc/'+bed_prefix+'/scores','class'+str(cell_class)+'_cluster'+str(cell_cluster)+'-'+phenotype_code+'.results')

log_file = re.sub('.results$','.log',results_file)

# Borrowed from https://github.com/shendurelab/mouse-atac/blob/master/ldscore_regression/gather_score_sumstats_results.py
log_entries = {}
for line in open(log_file):
	if 'SNPs with chi' in line:
		snp_count = int(re.search('\(([0-9]+) SNPs', line).group(1))
	if 'Total Observed scale h2' in line:
		match = re.search(': ([\-]*[0-9]\.[0-9]+(e[\-0-9]+)?) \((0\.[0-9]+)\)', line)
		h2 = float(match.group(1))
		h2_se = float(match.group(3))

# Combine meta data and log stats with partitioned h2 results
h2_results = pd.concat([
	pd.DataFrame(
		[[cell_cluster_name,phenotype_code,snp_count,h2,h2_se]],
		columns=['Cell_cluster','Phenotype','SNP_count','h2','h2_std_error']
	),
	pd.read_csv(results_file,header=0,index_col=0,sep='\t').loc[['L2_0']].reset_index()
],axis=1)

# Write in row name to make finding auxiliary files easier
h2_results.index = ['class'+str(cell_class)+'_cluster'+str(cell_cluster)+'-'+phenotype_code]

summary_file = re.sub('.results$','_summary.txt',results_file)

h2_results.to_csv(summary_file,sep='\t',header=True,index=True)

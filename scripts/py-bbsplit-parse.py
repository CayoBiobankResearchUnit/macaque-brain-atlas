#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import numpy as np
import gzip as gz
import os
import pyreadr as pyr
import pickle as pkl

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
# sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','NSM084','macaque','human','mouse']

prefix = arguments[1]
this_id = arguments[2]
genome1 = arguments[3]
genome2 = arguments[4]
genome3 = arguments[5]

adata = ad.read_h5ad('hdf5/anndata_rna-'+this_id+'.h5ad')

annotations = pd.read_csv('stats/clusters/rna-final-cellclasses.txt',sep='\t',header=None,index_col=None)

annotations.index = annotations[0].to_list()

annotations = annotations.rename(columns={1:'cell_class'})[['cell_class']]

meta = adata.obs[['id','cell','n.umi']]

meta = meta.merge(annotations,left_index=True,right_index=True,how='left')

# meta['n_'+genome1] = 0
# meta['n_'+genome2] = 0
# meta['n_ambiguous'] = 0

# <(samtools view bam_region/${this_id}_region_${chr}_${pos}.bam | cut -f 1 | cut -d '|' -f 3-5 | sed 's/|/_/g' | sed 's/^/'${this_id}'-/g') \

genome1_indexes = []
with open('fastq/bbsplit/'+this_id+'_'+genome1+'.fq', 'r') as f:
	for count, line in enumerate(f, start=0):
		if count % 4 == 0:
			# print(line)
			this_cell = this_id + '-'+'_'.join(line.split('|')[2:5])
			genome1_indexes.append(this_cell)

genome2_indexes = []
with open('fastq/bbsplit/'+this_id+'_'+genome2+'.fq', 'r') as f:
	for count, line in enumerate(f, start=0):
		if count % 4 == 0:
			# print(line)
			this_cell = this_id + '-'+'_'.join(line.split('|')[2:5])
			genome2_indexes.append(this_cell)

genome3_indexes = []
with open('fastq/bbsplit/'+this_id+'_'+genome3+'.fq', 'r') as f:
	for count, line in enumerate(f, start=0):
		if count % 4 == 0:
			# print(line)
			this_cell = this_id + '-'+'_'.join(line.split('|')[2:5])
			genome3_indexes.append(this_cell)


genome1_counts = pd.DataFrame(pd.DataFrame({0:genome1_indexes})[0].value_counts()).rename(columns={0:'n_'+genome1})
genome2_counts = pd.DataFrame(pd.DataFrame({0:genome2_indexes})[0].value_counts()).rename(columns={0:'n_'+genome2})
genome3_counts = pd.DataFrame(pd.DataFrame({0:genome3_indexes})[0].value_counts()).rename(columns={0:'n_'+genome3})

all_counts = genome1_counts.merge(genome2_counts,how='outer',left_index=True,right_index=True).merge(genome3_counts,how='outer',left_index=True,right_index=True)

all_counts.loc[[np.isnan(i) for i in all_counts['n_'+genome1]],'n_'+genome1] = 0
all_counts.loc[[np.isnan(i) for i in all_counts['n_'+genome2]],'n_'+genome2] = 0
all_counts.loc[[np.isnan(i) for i in all_counts['n_'+genome3]],'n_'+genome3] = 0

all_counts['n_'+genome1] = all_counts['n_'+genome1].astype(int)
all_counts['n_'+genome2] = all_counts['n_'+genome2].astype(int)
all_counts['n_'+genome3] = all_counts['n_'+genome3].astype(int)

all_counts['n_assigned'] = all_counts.sum(axis=1).to_list()

meta = meta.merge(all_counts,how='outer',left_index=True,right_index=True)

meta['id'] = this_id
meta['cell'] = [x.replace(this_id+'-','') for x in meta.index.to_list()]

meta['frac_'+genome1] = meta['n_'+genome1] / meta['n_assigned']
meta['frac_'+genome2] = meta['n_'+genome2] / meta['n_assigned']
meta['frac_'+genome3] = meta['n_'+genome3] / meta['n_assigned']

meta['rt'] = meta['cell'].str.split(r'[_]').map(lambda x: x[2]).to_list()

meta['rt'].value_counts()

meta.groupby('rt')['frac_'+genome2].median()
meta.groupby('rt')['frac_'+genome2].mean()

meta[(meta['frac_'+genome1] > 0.3333) & (meta['n_assigned'] >= 25)]['cell_class'].value_counts()
meta[(meta['frac_'+genome2] > 0.3333) & (meta['n_assigned'] >= 25)]['cell_class'].value_counts()
meta[(meta['frac_'+genome3] > 0.3333) & (meta['n_assigned'] >= 25)]['cell_class'].value_counts()

os.makedirs('stats/contamination/'+genome1+'_'+genome2+'_'+genome3,exist_ok=True)

meta.to_csv('stats/contamination/'+genome1+'_'+genome2+'_'+genome3+'/'+this_id+'_assigned.txt.gz',sep='\t',index=None,header=True)
meta.to_pickle('stats/contamination/'+genome1+'_'+genome2+'_'+genome3+'/'+this_id+'_assigned.pkl')

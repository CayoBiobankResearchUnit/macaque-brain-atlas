#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import numpy as np
import gzip as gz
import scrublet as scr
import os

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-prep.py','atac','atac','1']

prefix = arguments[1]
analysis = arguments[2]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_table('data/'+prefix+'_metadata.txt').id.tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[int(arguments[3])-1]

# anndata matrices are transposed relative to what we usually use in monocle3
mm = si.mmread('mm/'+prefix+'-'+query_sample+'-expr.mm').T.tocsc()                                                       n_prin_comps=30)

# Initialize a vector of cell names and read them in from a gzipped file (splitlines to remove trailing newline)
cells = []
for line in gz.open('mm/'+prefix+'-'+query_sample+'-expr.cols.txt.gz',mode='rt',encoding='utf8', newline='\n').read().splitlines():
	cells.append(line)

# Do the same thing with features
features = []
for line in gz.open('mm/'+prefix+'-'+query_sample+'-expr.rows.txt.gz',mode='rt',encoding='utf8', newline='\n').read().splitlines():
	features.append(line)

# Read in the data frames for the cell metadata and the feature metadata
meta_cells = pd.read_table('mm/'+prefix+'-'+query_sample+'-metadata.txt.gz',compression='gzip')
meta_features = pd.read_table('mm/'+prefix+'-'+query_sample+'-features.txt.gz',compression='gzip')
meta_samples = pd.read_table('data/'+prefix+'_run-data.txt.gz',compression='gzip')

meta_sample = meta_samples[meta_samples['sample']==query_sample]
collision_rate = meta_sample['collision_rate'].squeeze()
meta_cells['sequencing_run_id'] = meta_sample['run_id'].squeeze()

meta_cells.index = cells
meta_features.index = features

# Create an AnnData object from the sparse matrix
adata = ad.AnnData(mm,dtype='int',obs=meta_cells,var=meta_features)

# Write in all the relevant slots
adata.obs_names = pd.Index(cells)
adata.var_names = pd.Index(features)

# Filter
if analysis == 'atac':
	adata = adata[(adata.obs.umi_binarized>=1000) & (adata.obs.umi_binarized <= 1e5) & (adata.obs.FRIP>=0.3)]

# adata.obs = adata.obs.reset_index(drop=True)

try:
	adata = sc.external.pp.scrublet(adata,expected_doublet_rate=0.05,copy=True)
except:
	adata = sc.external.pp.scrublet(adata,expected_doublet_rate=0.05,copy=True,threshold=1.0)

# Make a Scrublet score distribution plot
sc.external.pl.scrublet_score_distribution(adata,save='_'+prefix+'-'+query_sample+'.pdf')

# By default the above writes the file "scrublet_score_distribution" with the suffix above appended
os.makedirs('figures/scrublet',exist_ok=True)
os.rename('figures/scrublet_score_distribution_'+prefix+'-'+query_sample+'.pdf','figures/scrublet/'+'scrublet_score_distribution_'+prefix+'-'+query_sample+'.pdf')

# Save Scrublet stats for custom plotting
os.makedirs('stats/scrublet',exist_ok=True)

# Write simulated scrublet scores
#with open('stats/scrublet/scrublet_score_distribution_'+prefix+'-'+query_sample+'-sim.txt',mode='wt',encoding='utf8') as f:
with gz.open('stats/scrublet/scrublet_score_distribution_'+prefix+'-'+query_sample+'-sim.txt.gz',mode='wt',encoding='utf8') as f:
	for line in adata.uns['scrublet']['doublet_scores_sim'].tolist():
		f.write(str(line))
		f.write('\n')

# Write observed scrublet scores
#with open('stats/scrublet/scrublet_score_distribution_'+prefix+'-'+query_sample+'-obs.txt',mode='wt',encoding='utf8') as f:
with gz.open('stats/scrublet/scrublet_score_distribution_'+prefix+'-'+query_sample+'-obs.txt.gz',mode='wt',encoding='utf8') as f:
	for line in adata.obs.doublet_score.tolist():
		f.write(str(line))
		f.write('\n')

# Write scrublet threshold
scrublet_parameters = {**{'id':query_sample}, **adata.uns['scrublet']['parameters'], **{'threshold':adata.uns['scrublet']['threshold']}}

pd.DataFrame(scrublet_parameters,index=[0]).to_csv('stats/scrublet/scrublet_parameters_'+prefix+'-'+query_sample+'.txt',sep='\t',index=False)

# Make an hdf5 directory if needed
os.makedirs('hdf5',exist_ok=True)

# adata.X = adata.X.astype(np.int32)

# Write the AnnData object to an hdf5 file
adata.write_h5ad(filename='hdf5/anndata_'+prefix+'-'+query_sample+'.h5ad')

# Write the AnnData object to an hdf5 file
# adata.write_loom(filename='hdf5/anndata_'+prefix+'-'+query_sample+'.loom')
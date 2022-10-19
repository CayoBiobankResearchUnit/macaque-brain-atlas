#!/usr/bin/env python
# 🐛

import sys
import scipy
import scipy.cluster.hierarchy
import anndata as ad
import scanpy as sc
import numpy as np
import os
import pandas as pd
import numpy as np

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered_endogenous.h5ad')

cellclasses = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
cellclass_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

cellclasses[1] = pd.Categorical(
	values=cellclasses[1],
	categories=cellclass_levels,
)

cellclasses = cellclasses.rename(columns={0:'cell',1:'cell_class'})

meta = adata.obs.merge(cellclasses,how='left',on='cell')

missing_cells = meta['cell_class'].value_counts()[meta['cell_class'].value_counts() == 0].index.to_list()

meta['cell_class'] = meta['cell_class'].cat.remove_categories(missing_cells)
meta.index = meta['cell'].to_list()
adata.obs = meta

cellsubclusters = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters.txt',header=None,sep='\t')
cellsubcluster_levels = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

cellsubclusters[1] = pd.Categorical(
	values=cellsubclusters[1],
	categories=cellsubcluster_levels,
)

cellsubclusters = cellsubclusters.rename(columns={0:'cell',1:'cell_subcluster'})

meta = adata.obs.merge(cellsubclusters,how='left',on='cell')

missing_cells = meta['cell_subcluster'].value_counts()[meta['cell_subcluster'].value_counts() == 0].index.to_list()

meta['cell_subcluster'] = meta['cell_subcluster'].cat.remove_categories(missing_cells)
meta.index = meta['cell'].to_list()
adata.obs = meta

adata.obs['region'] = pd.Categorical(
	adata.obs['region'].astype('str'),
	categories = pd.Index(adata.obs['region'].cat.categories.to_list() + ['MT']).sort_values().to_list()
)

adata.obs['region'][(adata.obs['region'] == 'STS') & (adata.obs['animal_id'] == '2C0')] = 'MT'

# Fix cell class
keep_cells = adata.obs['cell_class'].value_counts()[adata.obs['cell_class'].value_counts() > 2000].index.to_list()
keep_categories = adata.obs['cell_class'].cat.categories[[i in keep_cells for i in adata.obs['cell_class'].cat.categories]].to_list()

adata.obs['cell_class'] = pd.Categorical(
	adata.obs['cell_class'].astype('str'),
	categories = keep_categories
)

n_min = 1000

out_counts = []
for j in range(0,n_min):
	out_counts.append([])

# out_counts.update({i:[]})

for i in adata.obs['cell_class'].cat.categories.to_list():
	print(i)
	this = adata[adata.obs['cell_class'] == i]
	for j in range(0,n_min):
		print(j)
		this_sampled = this[np.random.choice(this.n_obs, size=n_min, replace=False)]
		out_counts[j].append(np.sum(this_sampled.X,axis=0).tolist()[0])

os.makedirs('expr_bootstrap',exist_ok=True)
for j in range(0,n_min):
	print(j)
	x = np.array(out_counts[j]).astype('int')
	np.savetxt('expr_bootstrap/'+prefix+'_cellclass_pseudobulk_'+str(j+1).zfill(4)+'.txt.gz',x,delimiter='\t',fmt='%i')

# Save cells and genes
pd.DataFrame(adata.obs['cell_class'].cat.categories.to_list()).to_csv('expr_bootstrap/'+prefix+'_cellclass_pseudobulk_cells.txt.gz',header=None,index=None)
pd.DataFrame(adata.var.index.to_list()).to_csv('expr_bootstrap/'+prefix+'_cellclass_pseudobulk_genes.txt.gz',header=None,index=None)

# One final pseudobulk without sampling
out_counts = []

for i in adata.obs['cell_class'].cat.categories.to_list():
	print(i)
	this = adata[adata.obs['cell_class'] == i]
	out_counts.append(np.sum(this.X,axis=0).tolist()[0])

x = np.array(out_counts).astype('int')
np.savetxt('expr_bootstrap/'+prefix+'_cellclass_pseudobulk_'+str(0).zfill(4)+'.txt.gz',x,delimiter='\t',fmt='%i')
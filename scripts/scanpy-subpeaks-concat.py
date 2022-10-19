#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import os
import scglue

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-subpeaks-concat.py','15']

prefix = 'atac'
analysis = 'atac'

cell_class = int(arguments[1])-1

celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

cell = celltype_levels[cell_class]

celltypes_class = celltypes[(celltypes['glue_type'] == cell) & (celltypes['glue_type_confidence'] >= 0.95)]

sample_list = list(celltypes_class['id'].unique())

all_h5ad = []

for i in sample_list:
	sys.stderr.write('Now appending '+str(i)+' for cell class '+str(cell_class+1)+' ('+cell+')\n')
	this = ad.read_h5ad('hdf5/subpeaks/anndata_'+prefix+'-class'+str(cell_class+1)+'-'+str(i)+'.h5ad')
	all_h5ad.append(this)

# for i in range(0,len(all_h5ad)):
# 	all_h5ad[i].X = all_h5ad[i].X.astype(np.int64)

adata = ad.concat(all_h5ad,axis=0,join='inner',merge='same')

# Calculate new UMI
# good_cells = list(np.where(list(np.array(np.sum(mm,axis=1)).reshape(-1,)))[0])

adata.obs['subpeak_total_umi'] = list(np.array(np.sum(adata.X,axis=1)).reshape(-1,))
adata.var['subpeak_umi_total'] = list(np.array(np.sum(adata.X,axis=0)).reshape(-1,))
adata.var['subpeak_umi_per_cell'] = adata.var['subpeak_umi_total'] / adata.n_obs

# Binarize matrix
adata.X = (adata.X > 0).astype('int')

adata.obs['subpeak_binarized_umi'] = list(np.array(np.sum(adata.X,axis=1)).reshape(-1,))
adata.var['subpeak_umi_binarized'] = list(np.array(np.sum(adata.X,axis=0)).reshape(-1,))
adata.var['subpeak_umi_binarized_per_cell'] = adata.var['subpeak_umi_binarized'] / adata.n_obs

# Do an LSI while we're at it and then it's preprocessed
scglue.data.lsi(adata, n_components=100, n_iter=15)

# Write with a new prefix to facilitate future analysis
at_prefix=prefix+'sub'
adata.write_h5ad(filename='hdf5/subintegration/anndata_'+at_prefix+'_class'+str(cell_class+1)+'_preprocessed.h5ad')

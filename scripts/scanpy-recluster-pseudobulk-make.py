#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import numpy as np
import gzip as gz
import os
import matplotlib.pyplot as plt
import math
import pyreadr as pyr

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','rna','rna','3']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
push_status(prefix+' read_h5ad cluster_'+str(this_cluster+1))

# Read in cell annotations
# celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellcluster2s.txt',header=None,sep='\t')
# celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellcluster1s.txt',header=None,sep='\t')
celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellclusters.txt',header=None,sep='\t')
celltype_levels = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

celltype_levels = [i for i in celltype_levels if i in celltypes[1].to_list()]

adata.obs['cell_subcluster'] = pd.Categorical(
	values=celltypes[1],
	categories = celltype_levels
)

# Pseudobulk (at cell_subcluster level, can also do cell_type, etc.)

n_min = adata.obs['cell_subcluster'].value_counts().min()

cell_types = celltype_levels

out_counts = []
out_cells = []
out_classes = []
out_umi = []
for i in cell_types:
	this = adata[adata.obs['cell_subcluster'] == i]
	this_sampled = this[np.random.choice(this.n_obs, size=n_min, replace=False)]
	out_cells.append(len(this))
	out_counts.append(np.sum(this.X,axis=0).tolist()[0])
	out_umi.append(np.sum(this.X,axis=None))

bdata = ad.AnnData(
	np.array(out_counts),
	dtype='int',
	obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells,'umi':out_umi,'n_sampled':n_min},index=cell_types),
	var = adata.var
)

os.makedirs('stats/pseudobulk',exist_ok=True)
pd.DataFrame(bdata.X.transpose()).set_axis(cell_types,axis=1).set_axis(bdata.var.index.to_list(),axis=0).to_csv('stats/pseudobulk/'+prefix+'_class'+str(this_cluster+1)+'_pseudobulk_min_extr.txt.gz',sep='\t',index=True,header=True)
bdata.obs.to_csv('stats/pseudobulk/'+prefix+'_class'+str(this_cluster+1)+'_pseudobulk_min_metadata.txt.gz',sep='\t',index=True,header=True)

n_min = 1000

cell_types = adata.obs['cell_subcluster'].cat.categories[adata.obs['cell_subcluster'].value_counts() >= 1000].to_list()

out_counts = []
out_cells = []
out_classes = []
out_umi = []
for i in cell_types:
	this = adata[adata.obs['cell_subcluster'] == i]
	this_sampled = this[np.random.choice(this.n_obs, size=n_min, replace=False)]
	out_cells.append(len(this))
	out_counts.append(np.sum(this.X,axis=0).tolist()[0])
	out_umi.append(np.sum(this.X,axis=None))

bdata = ad.AnnData(
	np.array(out_counts),
	dtype='int',
	obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells,'umi':out_umi,'n_sampled':n_min},index=cell_types),
	var = adata.var
)

pd.DataFrame(bdata.X.transpose()).set_axis(cell_types,axis=1).set_axis(bdata.var.index.to_list(),axis=0).to_csv('stats/pseudobulk/'+prefix+'_class'+str(this_cluster+1)+'_pseudobulk_1k_extr.txt.gz',sep='\t',index=True,header=True)
bdata.obs.to_csv('stats/pseudobulk/'+prefix+'_class'+str(this_cluster+1)+'_pseudobulk_1k_metadata.txt.gz',sep='\t',index=True,header=True)

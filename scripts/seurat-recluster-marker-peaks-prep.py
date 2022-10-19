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
import random

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-recluster-marker-peaks-prep.py','atacsub','subpeak','12']

prefix = arguments[1]
# analysis = arguments[2]
glue_prefix = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/subintegration/anndata_'+prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')

# Read in cell annotations
celltypes = pd.read_csv('stats/clusters/subintegration/'+glue_prefix+'-'+prefix+'-class'+str(this_cluster+1)+'-cellsubtype-predictions.txt',header=0,sep='\t',index_col=0)

celltypes = celltypes[['glue_subtype','glue_subtype_confidence']]

adata.obs = adata.obs.merge(celltypes,how='left',left_on='cell',right_index=True,copy=True)

# adata = adata[adata.obs['glue_subtype_confidence'] >= 0.5]

cell_subtypes = pd.read_csv('stats/subclusters/'+'rna-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

adata.obs['cell_subcluster'] = pd.Categorical(
	values = adata.obs['glue_subtype'],
	categories = cell_subtypes
)

adata.obs['cell_subcluster'] = adata.obs['cell_subcluster'].cat.remove_unused_categories()


# Cell sample threshold (keep only this many cells for each cell subtype)
cell_sample_threshold = 1000

cell_counts = adata.obs['cell_subcluster'].value_counts(sort=False)
cell_counts_to_sample = cell_counts.copy()

cell_counts_to_sample[cell_counts_to_sample > cell_sample_threshold] = cell_sample_threshold

keep_cells = []
for i in adata.obs['cell_subcluster'].cat.categories:
	if cell_counts[i] > cell_sample_threshold:
		print('Sampling '+i)
		n_to_sample = cell_counts_to_sample[i]
		cells = adata.obs[adata.obs['cell_subcluster']==i].index.to_list()
		random.seed(a = 1)
		cells_to_keep = random.sample(cells,cell_sample_threshold)
	else:
		cells_to_keep = adata.obs[adata.obs['cell_subcluster']==i].index.to_list()
	keep_cells = keep_cells + cells_to_keep

adata_out = adata[keep_cells]

os.makedirs('mm/recluster/seurat',exist_ok=True)
si.mmwrite(target='mm/recluster/seurat/'+prefix+'_downsampled_class'+str(this_cluster+1)+'_counts',a=adata_out.X)
adata_out.obs.to_csv('mm/recluster/seurat/'+prefix+'_downsampled_class'+str(this_cluster+1)+'_cell_metadata.txt.gz',sep='\t',index=True,header=True)
adata_out.var.to_csv('mm/recluster/seurat/'+prefix+'_downsampled_class'+str(this_cluster+1)+'_feature_metadata.txt.gz',sep='\t',index=True,header=True)

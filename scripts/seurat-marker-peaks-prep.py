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
# arguments = ['seurat-marker-peaks-prep.py','atac','atac']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered.h5ad')
push_status(prefix+' read_h5ad')

# Read in cell annotations
celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

celltypes = celltypes[['modality','glue_type','glue_type_confidence']]

celltypes = celltypes[celltypes['glue_type_confidence'] >= 0.95]

adata.obs = adata.obs.merge(celltypes,how='left',left_on='cell',right_index=True,copy=True)

adata.obs['glue_type'] = pd.Categorical(
	values = adata.obs['glue_type'],
	categories = celltype_levels
)

# adata = adata[adata.obs['glue_type_confidence'] > 0.95]

adata.obs['glue_type'] = adata.obs['glue_type'].cat.remove_unused_categories()
adata.obs['cell_class'] = adata.obs['glue_type']

cell_types = adata.obs['cell_class'].cat.categories.to_list()

adata = adata[~adata.obs['cell_class'].isnull()]

# Cell sample threshold (keep only this many cells for each cell subtype)
cell_sample_threshold = 1000

cell_counts = adata.obs['cell_class'].value_counts(sort=False)
cell_counts_to_sample = cell_counts.copy()

cell_counts_to_sample[cell_counts_to_sample > cell_sample_threshold] = cell_sample_threshold

keep_cells = []
for i in adata.obs['cell_class'].cat.categories:
	if cell_counts[i] > cell_sample_threshold:
		print('Sampling '+i)
		n_to_sample = cell_counts_to_sample[i]
		cells = adata.obs[adata.obs['cell_class']==i].index.to_list()
		random.seed(a = 1)
		cells_to_keep = random.sample(cells,cell_sample_threshold)
	else:
		cells_to_keep = adata.obs[adata.obs['cell_class']==i].index.to_list()
	keep_cells = keep_cells + cells_to_keep

adata_out = adata[keep_cells]

os.makedirs('mm/seurat',exist_ok=True)
si.mmwrite(target='mm/seurat/'+prefix+'_downsampled_counts',a=adata_out.X)
adata_out.obs.to_csv('mm/seurat/'+prefix+'_downsampled_cell_metadata.txt.gz',sep='\t',index=True,header=True)
adata_out.var.to_csv('mm/seurat/'+prefix+'_downsampled_feature_metadata.txt.gz',sep='\t',index=True,header=True)

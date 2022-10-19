#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import os
import matplotlib.pyplot as plt

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-prep.py','atac','atac']

prefix = arguments[1]
analysis = arguments[2]

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_geneactivity.h5ad')
push_status(prefix+' read_h5ad')

adata_final = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_postprocessed.h5ad')

# Temporary fix
adata.obs['hemisphere'] = pd.Categorical(adata.obs['hemisphere'], categories=['','L','R'], ordered=False)
adata.obs.iloc()[(adata.obs['region']=='CV').to_list(),(adata.obs.columns=='hemisphere')] = ''
adata.obs.iloc()[(adata.obs['region']=='MB').to_list(),(adata.obs.columns=='hemisphere')] = ''

# Count normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
push_status(prefix+' normalize_total')

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p')

# Add in the UMAP
adata.obsm['X_umap'] = adata_final.obsm['X_umap'].copy()

# Temporary
adata.obsm['X_umap'] = np.loadtxt('stats/umap_post/umap02_atac_all_mindist_0.25.txt.gz',dtype='float32')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
push_status(prefix+' write_h5ad')

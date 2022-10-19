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
import scrublet as scr
import os
import math

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-rna-postprocess.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# After detecting exogenous contamination, this script is to reprocess after removing them. It's for RNA only

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered.h5ad')
push_status(prefix+' read_h5ad')

# Read in classes
cellclasses = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
cellclass_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

cellclasses[1] = pd.Categorical(
	values=cellclasses[1],
	categories=cellclass_levels,
)

cellclasses = cellclasses.rename(columns={0:'cell',1:'cell_class'})

meta = adata.obs.merge(cellclasses,how='left',on='cell')
meta.index = meta['cell'].to_list()
adata.obs = meta

adata = adata[(adata.obs['cell_class'] != 'radial glial cells') & (adata.obs['cell_class'] != 'mesenchymal stem cells')]

adata.obs = adata.obs.drop(['cell_class'],axis=1)
push_status(prefix+' prepped')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_filtered_endogenous.h5ad')

# Rasterize plots
sc._settings.settings._vector_friendly=True

adata.var['mt'] = adata.var.chromosome == 'MT'

# Calculate QC metrics (including mitochondrial percentage)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
push_status(prefix+' calculate_qc_metrics')	
# Remove sex chromosomes, MT, and unplaced scaffolds
adata = adata[:,np.isin(adata.var['chromosome'],[str(i) for i in range(1,21,1)])]

# Save a copy of the unnormalized matrix (save to disk in order to conserve memory)
# adata.layers['raw'] = adata.X.copy()
os.makedirs('npz',exist_ok=True)
sp.sparse.save_npz('npz/'+prefix+'_all_raw_final_endogenous.npz',adata.X,compressed=True)
push_status(prefix+' save_npz raw')

# adata.X = adata.X.asfptype()

# Count normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
push_status(prefix+' normalize_total')

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p')

# Save a copy of the normalized matrix (save to disk in order to conserve memory)
# adata.layers['normalized'] = adata.X.copy()
sp.sparse.save_npz('npz/'+prefix+'_all_normalized_final_endogenous.npz',adata.X,compressed=True)
push_status(prefix+' save_npz normalized')

# Identify variable genes
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata, n_top_genes=10000)
push_status(prefix+' highly_variable_genes')

# Write var to file
os.makedirs('stats/var_final',exist_ok=True)
adata.var.to_csv('stats/var_final/var_'+prefix+'_all_endogenous.txt.gz',sep='\t',header=True,index=True)

# Subset to variable genes
adata = adata[:,adata.var['highly_variable']]
# sc.pp.highly_variable_genes(adata, n_top_genes=10000,subset=True)
push_status(prefix+' highly_variable_genes subset')

# Regress out technical variables
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.regress_out(adata, ['total_counts'])
push_status(prefix+' regress_out')

# Scale each gene
sc.pp.scale(adata, zero_center=True)
push_status(prefix+' scale')

# Do PCA (use highly variable genes)
sc.pp.pca(adata, n_comps=n_pcs, svd_solver='arpack', use_highly_variable=True)
push_status(prefix+' pca')

# Set number of downstream dimensions to n_pcs
n_dim = n_pcs

os.makedirs('stats/preintg_post',exist_ok=True)

np.savetxt('stats/preintg_post/pca_'+prefix+'_all_endogenous.txt.gz',adata.obsm['X_pca'],delimiter='\t')
push_status(prefix+' savetxt preintg')

# Write metadata
os.makedirs('stats/umap_post',exist_ok=True)
adata.obs.to_csv('stats/umap_post/meta_'+prefix+'_all_endogenous.txt.gz',sep='\t',header=True,index=True)

# Integrate with bbknn
sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.floor(15/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
push_status(prefix+' bbknn 2D')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_post/umap02_'+prefix+'_all_endogenous.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D')

# Redo neighbors and UMAP to make 10-dimensional UMAP
sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.ceil(30/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
push_status(prefix+' bbknn 10D')
# sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='euclidean')
# push_status(prefix+' neighbors 10D')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_post/umap10_'+prefix+'_all_endogenous.txt.gz',umap10,delimiter='\t')
push_status(prefix+' umap 10D')

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_endogenous.h5ad')
push_status(prefix+' write_h5ad final')

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_dim, method='umap', metric='cosine')
push_status(prefix+' neighbors 2D')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_post/umap02_'+prefix+'_all_endogenous_uncorrected.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D')

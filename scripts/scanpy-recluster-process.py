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
import episcanpy.api as epi

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','rna','rna','11']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
push_status(prefix+' read_h5ad'+' cluster_'+str(this_cluster+1))

if analysis == 'rna':
	
	adata.var['mt'] = adata.var.chromosome == 'MT'
	
	# Calculate QC metrics (including mitochondrial percentage)
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	push_status(prefix+' calculate_qc_metrics'+' cluster_'+str(this_cluster+1))

	adata = adata[:,np.isin(adata.var['chromosome'],[str(i) for i in range(1,21,1)])]
	
	# Count normalize the data
	sc.pp.normalize_total(adata, target_sum=1e4)
	push_status(prefix+' normalize_total'+' cluster_'+str(this_cluster+1))
	
	# Logarithmize the data
	sc.pp.log1p(adata)
	push_status(prefix+' log1p'+' cluster_'+str(this_cluster+1))
	
	# Identify variable genes
	# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	
	try:
		sc.pp.highly_variable_genes(adata, n_top_genes=2000)
	except IndexError:
		print('sc.pp.highly_variable_genes() returned an error. Proceed with n_top_genes=None.')
		sc.pp.highly_variable_genes(adata)
	
	push_status(prefix+' highly_variable_genes'+' cluster_'+str(this_cluster+1))
	
	# # sc.pp.highly_variable_genes(adata, n_top_genes=10000,subset=True)
	# push_status(prefix+' highly_variable_genes subset'+' cluster_'+str(this_cluster+1))
	
	# Subset to variable genes
	adata = adata[:,adata.var['highly_variable']]
	
	bdata_list = []
	# Split by batch
	for i in adata.obs['sequencing_run_id'].cat.categories:
		this = adata[adata.obs['sequencing_run_id'] == i].copy()
		sc.pp.regress_out(this, ['total_counts'])
		# sc.pp.scale(this, zero_center=True)
		bdata_list.append(this)
	
	bdata = ad.concat(bdata_list,axis=0,join='inner',merge='same',uns_merge='same')
	
	sc.pp.scale(bdata, zero_center=True)
	push_status(prefix+' scale'+' cluster_'+str(this_cluster+1))
	
	adata = bdata[adata.obs.index].copy()
	
	# Regress out technical variables
	# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
	
	# sc.pp.regress_out(adata, ['total_counts'])
	# push_status(prefix+' regress_out'+' cluster_'+str(this_cluster+1))
	# 
	# # Scale each gene
	# sc.pp.scale(adata, zero_center=True)
	# push_status(prefix+' scale'+' cluster_'+str(this_cluster+1))
	
	# Do PCA (use highly variable genes)
	if min(adata.shape) < (n_pcs - 1):
		n_pcs = min(adata.shape) - 1
	
	# sc.pp.pca(adata, n_comps=n_pcs, svd_solver='arpack', use_highly_variable=False)
	sc.pp.pca(adata, n_comps=n_pcs, svd_solver='arpack')
	push_status(prefix+' pca'+' cluster_'+str(this_cluster+1))
	
	# Set number of downstream dimensions to n_pcs
	n_dim = n_pcs

os.makedirs('stats/umap_recluster',exist_ok=True)

adata.obs.to_csv('stats/umap_recluster/meta_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)

try:
	sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.floor(15/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
except ValueError:
	print('sce.pp.bbknn() returned an error. Proceed with sc.pp.neighbors().')
	sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_dim, method='umap', metric='cosine')
	push_status(prefix+' neighbors 2D'+' cluster_'+str(this_cluster+1))
else:
	push_status(prefix+' bbknn 2D'+' cluster_'+str(this_cluster+1))

# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_dim, method='umap', metric='cosine')
# push_status(prefix+' neighbors 2D'+' cluster_'+str(this_cluster+1))

np.savetxt('stats/umap_recluster/pca_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',adata.obsm['X_pca'],delimiter='\t')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_recluster/umap02_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D'+' cluster_'+str(this_cluster+1))

# Redo neighbors and UMAP to make 10-dimensional UMAP
try:
	sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.ceil(30/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
except ValueError:
	print('sce.pp.bbknn() returned an error. Proceed with sc.pp.neighbors().')
	sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='cosine')
	push_status(prefix+' neighbors 10D'+' cluster_'+str(this_cluster+1))
else:
	push_status(prefix+' bbknn 10D'+' cluster_'+str(this_cluster+1))

# sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='cosine')
# push_status(prefix+' neighbors 10D'+' cluster_'+str(this_cluster+1))

os.makedirs('stats/connectivities',exist_ok=True)
np.save('stats/connectivities/connectivities_'+prefix+'_class'+str(this_cluster+1)+'.npy', adata.obsp['connectivities'])
np.savetxt('stats/connectivities/cellids_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',adata.obs.index.to_numpy(),delimiter='\n',fmt='%s')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_recluster/umap10_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',umap10,delimiter='\t')
push_status(prefix+' umap 10D'+' cluster_'+str(this_cluster+1))

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'_postprocessed.h5ad')
push_status(prefix+' write_h5ad final'+' cluster_'+str(this_cluster+1))

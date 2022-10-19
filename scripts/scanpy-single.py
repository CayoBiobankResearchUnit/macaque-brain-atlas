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
import episcanpy.api as epi

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-prep.py','rna','rna','1']

prefix = arguments[1]
analysis = arguments[2]

# Read in updated scrublet scores
scrublet_scores = pd.read_csv('data/'+prefix+'_scrublet_parameters_final.txt',sep='\t')

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[int(arguments[3])-1]

adata = ad.read_h5ad('hdf5/anndata_'+prefix+'-'+query_sample+'.h5ad')

adata.obs.insert(adata.obs.columns.get_loc('predicted_doublet')+1,'manual_doublet',(adata.obs.doublet_score > scrublet_scores[scrublet_scores['id']==query_sample]['final_threshold'].squeeze()).tolist()) 

adata = adata[~(adata.obs['manual_doublet'])]

if analysis == 'atac':
	
	epi.tl.find_genes(adata,
		gtf_file=gtf_file,
		key_added='transcript_annotation',
		upstream=2000,
		feature_type='transcript',
		annotation='ensembl',
		raw=False)
	
	# Remove any potential empty features or barcodes
	epi.pp.filter_cells(adata, min_features=1)
	epi.pp.filter_features(adata, min_cells=1)
	
	# Quality controls
	adata.obs['log_nb_features'] = [np.log10(x) for x in adata.obs['nb_features']]
	
	# Filter cells
	epi.pp.filter_cells(adata, min_features=min_features)
	epi.pp.filter_features(adata, min_cells=min_cells)
	
	# Calculate variable features
	epi.pp.cal_var(adata)
	
	# Do all dim reductions
	epi.pp.lazy(adata)
	
	umap_pre = adata.obsm['X_umap']
	pca_pre = adata.obsm['X_pca']
	tsne_pre = adata.obsm['X_tsne']
	
	adata.layers['binary'] = adata.X.copy()
	
	# The code below launches an error, but normalizes anyway
	try:
		epi.pp.normalize_total(adata)
	except:
		if adata.X.max() > 1:
			print('epi.pp.normalize_total() returned an error but normalized correctly. Proceed.')
		else:
			sys.exit('epi.pp.normalize_total() failed to normalize.')
	
	adata.layers['normalized'] = adata.X.copy()
	
	# Redo PCA, t-SNE, and UMAP
	epi.pp.lazy(adata,
		n_neighbors=15,
		min_dist=0.5,
		n_components=2,
		copy=False)
	
	# For clustering, save a copy
	bdata = epi.pp.lazy(adata,
		n_neighbors=30,
		min_dist=0,
		n_components=10,
		copy=True)
	
	# Save the 10-dimensional UMAP
	adata.obsm['X_umap10'] = bdata.obsm['X_umap']
	
elif analysis == 'rna':
	
	adata.var['mt'] = adata.var.chromosome == 'MT'
	
	# Calculate QC metrics (including mitochondrial percentage)
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	
	# Filter adata
	adata = adata[(adata.obs.total_counts >= 100) & (adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :]
	
	# Count normalize the data
	sc.pp.normalize_total(adata, target_sum=1e4)
	
	# Logarithmize the data
	sc.pp.log1p(adata)
	
	# Identify variable genes
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	
	# Save a copy of the raw matrix
	adata.raw = adata
	
	# Regress out technical variables
	sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
	
	# Scale each gene
	sc.pp.scale(adata, max_value=10)
	
	# Do PCA (use highly variable)
	sc.pp.pca(adata, n_comps=50, svd_solver='arpack', use_highly_variable=True)
	
	# Compute neighborhood graph
	sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, method='umap', metric='euclidean')
	
	# Do UMAP
	sc.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2)
	
	# Make a copy
	bdata = adata.copy()
	
	# Redo neighborhood graph and UMAP with adjusted parameters for clustering
	sc.pp.neighbors(bdata, n_neighbors=30, n_pcs=50, method='umap', metric='euclidean')
	sc.tl.umap(bdata, min_dist=0, spread=1.0, n_components=10)
	
	# Save the 10-dimensional UMAP
	adata.obsm['X_umap10'] = bdata.obsm['X_umap']

# Make sure UMAP folder exists
os.makedirs('stats/umap',exist_ok=True)

# Write the 2D dimensional and 10D UMAP coordinates to a file

adata.obs.to_csv('stats/umap/meta_'+prefix+'-'+query_sample+'.txt.gz',sep='\t',header=True,index=True)
np.savetxt('stats/umap/umap02_'+prefix+'-'+query_sample+'.txt.gz',adata.obsm['X_umap'],delimiter='\t')
np.savetxt('stats/umap/umap10_'+prefix+'-'+query_sample+'.txt.gz',adata.obsm['X_umap10'],delimiter='\t')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'-'+query_sample+'_preprocessed.h5ad')
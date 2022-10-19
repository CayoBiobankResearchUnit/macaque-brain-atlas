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
import sklearn as sk
import sklearn.decomposition # TruncatedSVD
import sklearn.preprocessing # normalize
import sklearn.feature_extraction.text # TfidfTransformer

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

singlet_cells = pd.read_csv('stats/cells/'+prefix+'-singlets.txt',header=None,sep='\t')

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all.h5ad')
push_status(prefix+' read_h5ad')

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

# Temporary fix (midline structures were incorrectly assigned to the right hemisphere)
adata.obs['hemisphere'] = pd.Categorical(adata.obs['hemisphere'], categories=['','L','R'], ordered=False)
adata.obs.iloc()[(adata.obs['region']=='CV').to_list(),(adata.obs.columns=='hemisphere')] = ''
adata.obs.iloc()[(adata.obs['region']=='MB').to_list(),(adata.obs.columns=='hemisphere')] = ''

adata.obs['is_singlet'] = np.isin(adata.obs['cell'].to_list(),singlet_cells[0].to_list()).tolist()
adata = adata[adata.obs['is_singlet']]
push_status(prefix+' remove_doublet')

# Rasterize plots
sc._settings.settings._vector_friendly=True

if analysis == 'atac':
	
	# Annotate genes using GTF file
	epi.tl.find_genes(adata,
		gtf_file=gtf_file,
		key_added='transcript_annotation',
		upstream=2000,
		feature_type='transcript',
		annotation='ensembl',
		raw=False)
	push_status(prefix+' find_genes')
	
	# No need to filter cells
	# epi.pp.filter_cells(adata, min_features=min_features)
	# push_status(prefix+' filter_cells')
	
	epi.pp.filter_features(adata, min_cells=min_cells)
	push_status(prefix+' filter_features')
	
	# Remove sex chromosomes
	peaks_annotation = np.array([line.split('_') for line in adata.var_names.tolist()])
	adata.var['chromosome'] = peaks_annotation[:,0]
	adata.var['bp1'] = peaks_annotation[:,1]
	adata.var['bp2'] = peaks_annotation[:,2]
	adata = adata[:,np.isin(adata.var['chromosome'],[str(i) for i in range(1,21,1)])]
	
	# In case removal of sex chromosomes results in zero features, remove those cells
	epi.pp.filter_cells(adata, min_features=0)
	push_status(prefix+' filter_cells')
	
	# Quality controls
	adata.obs['log_nb_features'] = [np.log10(x) for x in adata.obs['nb_features']]
	
	# Calculate variable features
	epi.pp.cal_var(adata)
	push_status(prefix+' cal_var')
	
	# Annotate highly variable features
	adata.var['highly_variable'] = adata.var['variability_score'] >= adata.var['variability_score'].sort_values(ascending=False)[range(0,100000)].min()
	push_status(prefix+' select_var_feature')
	
	# Write var to file
	os.makedirs('stats/var_final',exist_ok=True)
	adata.var.to_csv('stats/var_final/var_'+prefix+'_all.txt.gz',sep='\t',header=True,index=True)
	
	# Get error if not floating-point type
	# https://stackoverflow.com/questions/8650014/sparse-matrix-valueerror-matrix-type-must-be-f-d-f-or-d
	#
	adata.X = adata.X.asfptype()
	
	# TF-IDF normalization (on highly variable)
	tfidf = sk.feature_extraction.text.TfidfTransformer(norm='l1')
	tfidf.fit(adata.X[:,adata.var['highly_variable']])
	tfidf_matrix = tfidf.transform(adata.X[:,adata.var['highly_variable']])
	push_status(prefix+' tfidf raw')
	
	# LSI
	svd_model = sk.decomposition.TruncatedSVD(n_components=n_pcs,algorithm='arpack',n_iter=10,random_state=42)
	adata.obsm['X_lsi'] = svd_model.fit_transform(tfidf_matrix)
	push_status(prefix+' lsi raw')
	
	# L2 normalization of PCA
	adata.obsm['X_pca'] = sk.preprocessing.normalize(adata.obsm['X_lsi'][:,range(1,n_pcs)],norm='l2')
	
	# Set number of downstream dimensions to n_pcs - 1 (dropping first PC)
	n_dim = n_pcs - 1
elif analysis == 'rna':
	
	adata.var['mt'] = adata.var.chromosome == 'MT'
	
	# Calculate QC metrics (including mitochondrial percentage)
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	push_status(prefix+' calculate_qc_metrics')
	
	# # Filter adata
	# adata = adata[(adata.obs.total_counts >= 100) & (adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :]
	
	# Remove sex chromosomes, MT, and unplaced scaffolds
	adata = adata[:,np.isin(adata.var['chromosome'],[str(i) for i in range(1,21,1)])]
	
	# Save a copy of the unnormalized matrix (save to disk in order to conserve memory)
	# adata.layers['raw'] = adata.X.copy()
	os.makedirs('npz',exist_ok=True)
	sp.sparse.save_npz('npz/'+prefix+'_all_raw_final.npz',adata.X,compressed=True)
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
	sp.sparse.save_npz('npz/'+prefix+'_all_normalized_final.npz',adata.X,compressed=True)
	push_status(prefix+' save_npz normalized')
	
	# Identify variable genes
	# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	sc.pp.highly_variable_genes(adata, n_top_genes=10000)
	push_status(prefix+' highly_variable_genes')
	
	# Write var to file
	os.makedirs('stats/var_final',exist_ok=True)
	adata.var.to_csv('stats/var_final/var_'+prefix+'_all.txt.gz',sep='\t',header=True,index=True)
	
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

np.savetxt('stats/preintg_post/pca_'+prefix+'_all.txt.gz',adata.obsm['X_pca'],delimiter='\t')
push_status(prefix+' savetxt preintg')

# Write metadata
os.makedirs('stats/umap_post',exist_ok=True)
adata.obs.to_csv('stats/umap_post/meta_'+prefix+'_all.txt.gz',sep='\t',header=True,index=True)

# Integrate with bbknn
sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.floor(15/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
push_status(prefix+' bbknn 2D')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_post/umap02_'+prefix+'_all.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D')

# Redo neighbors and UMAP to make 10-dimensional UMAP
sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.ceil(30/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
push_status(prefix+' bbknn 10D')
# sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='euclidean')
# push_status(prefix+' neighbors 10D')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_post/umap10_'+prefix+'_all.txt.gz',umap10,delimiter='\t')
push_status(prefix+' umap 10D')

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_postprocessed.h5ad')
push_status(prefix+' write_h5ad final')

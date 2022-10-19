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
import episcanpy as epi
import sklearn as sk
import sklearn.decomposition # TruncatedSVD
import sklearn.preprocessing # normalize
import sklearn.feature_extraction.text # TfidfTransformer

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','subpeak','atacsub','6']

prefix = arguments[1]
at_prefix = arguments[2]
cell_class = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/subintegration/anndata_'+at_prefix+'_class'+str(cell_class+1)+'_preprocessed.h5ad')
push_status(prefix+' read_h5ad cluster_'+str(cell_class+1))

# Bring in GLUE subtype predictions
cell_subclusters = pd.read_csv('stats/clusters/subintegration/'+prefix+'-'+at_prefix+'-class'+str(cell_class+1)+'-cellsubtype-predictions.txt',header=0,index_col=0,sep='\t')
new_obs = adata.obs.merge(cell_subclusters[['glue_subtype','glue_subtype_confidence']],left_index=True,right_index=True,how='left')
adata.obs = new_obs

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

# # Write var to file
# os.makedirs('stats/var_final',exist_ok=True)
# adata.var.to_csv('stats/var_final/var_'+prefix+'_all.txt.gz',sep='\t',header=True,index=True)

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


os.makedirs('stats/umap_recluster',exist_ok=True)

adata.obs.to_csv('stats/umap_recluster/meta_'+prefix+'_class'+str(cell_class+1)+'.txt.gz',sep='\t',header=True,index=True)

try:
	sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.floor(15/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
except ValueError:
	print('sce.pp.bbknn() returned an error. Proceed with sc.pp.neighbors().')
	sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_dim, method='umap', metric='cosine')
	push_status(prefix+' neighbors 2D'+' cluster_'+str(cell_class+1))
else:
	push_status(prefix+' bbknn 2D'+' cluster_'+str(cell_class+1))

# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_dim, method='umap', metric='cosine')
# push_status(prefix+' neighbors 2D'+' cluster_'+str(cell_class+1))

np.savetxt('stats/umap_recluster/pca_'+prefix+'_class'+str(cell_class+1)+'.txt.gz',adata.obsm['X_pca'],delimiter='\t')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_recluster/umap02_'+prefix+'_class'+str(cell_class+1)+'.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D'+' cluster_'+str(cell_class+1))

# Redo neighbors and UMAP to make 10-dimensional UMAP
try:
	sce.pp.bbknn(adata,batch_key='sequencing_run_id',use_rep='X_pca',neighbors_within_batch=math.ceil(30/len(adata.obs['sequencing_run_id'].unique())),approx=True,n_pcs=n_dim,use_annoy=False,metric='cosine')
except ValueError:
	print('sce.pp.bbknn() returned an error. Proceed with sc.pp.neighbors().')
	sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='cosine')
	push_status(prefix+' neighbors 10D'+' cluster_'+str(cell_class+1))
else:
	push_status(prefix+' bbknn 10D'+' cluster_'+str(cell_class+1))

# sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_dim, method='umap', metric='cosine')
# push_status(prefix+' neighbors 10D'+' cluster_'+str(cell_class+1))

os.makedirs('stats/connectivities',exist_ok=True)
np.save('stats/connectivities/connectivities_'+prefix+'_class'+str(cell_class+1)+'.npy', adata.obsp['connectivities'])
np.savetxt('stats/connectivities/cellids_'+prefix+'_class'+str(cell_class+1)+'.txt.gz',adata.obs.index.to_numpy(),delimiter='\n',fmt='%s')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/umap_recluster/umap10_'+prefix+'_class'+str(cell_class+1)+'.txt.gz',umap10,delimiter='\t')
push_status(prefix+' umap 10D'+' cluster_'+str(cell_class+1))

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(cell_class+1)+'_postprocessed.h5ad')
push_status(prefix+' write_h5ad final'+' cluster_'+str(cell_class+1))

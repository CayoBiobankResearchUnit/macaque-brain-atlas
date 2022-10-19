#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import harmonypy as hp
import numpy as np
import gzip as gz
import scrublet as scr
import os
import episcanpy.api as epi

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-prep.py','atac','atac','1']

prefix = arguments[1]
n_iter = arguments[3]

# adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_prenorm.h5ad')

# pca_pre = pd.read_csv('stats/prenorm/prenorm_pca_'+prefix+'_all.txt.gz',header=None,sep='\t')
umap_prenorm = pd.read_csv('stats/prenorm/prenorm_umap_'+prefix+'_all.txt.gz',header=None,sep='\t').to_numpy()

# Calculate UMAP using the previous UMAP to anchor
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_dimreduced.h5ad')
push_status('read_h5ad dimreduced n'+str(n_iter))

adata.obsm['X_pca'] = np.loadtxt('stats/preintg/pca_'+prefix+'_all.txt.gz',dtype='float32')

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, method='umap', metric='euclidean')
push_status('neighbors dimreduced n'+str(n_iter))

sc.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2, init_pos=umap_prenorm)
push_status('umap dimreduced n'+str(n_iter))

umap_preharm = adata.obsm['X_umap'].copy()

os.makedirs('stats/harmony',exist_ok=True)

# These files below in theory should be the same between iterations (confirm by checking md5 checksums of uncompressed files)
np.savetxt('stats/harmony/umap_prenorm_'+prefix+'_n'+str(n_iter)+'_all.txt.gz',umap_prenorm,delimiter='\t')
np.savetxt('stats/harmony/umap_preharm_'+prefix+'_n'+str(n_iter)+'_all.txt.gz',umap_preharm,delimiter='\t')

# Now do harmony

sce.pp.harmony_integrate(adata,key='sequencing_run_id',basis='X_pca',adjusted_basis='X_pca_harmony',max_iter_harmony=int(n_iter))
push_status('harmony_integrate n'+str(n_iter))

# Write the harmony PCA into the main slot
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony'].copy()

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, method='umap', metric='euclidean')
push_status('neighbors harmony n'+str(n_iter))

sc.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2, init_pos=umap_preharm)
push_status('umap harmony n'+str(n_iter))

umap_postharm = adata.obsm['X_umap'].copy()

np.savetxt('stats/harmony/umap_postharm_'+prefix+'_n'+str(n_iter)+'_all.txt.gz',umap_postharm,delimiter='\t')

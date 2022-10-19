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
import networkx as nx
import scglue
from matplotlib import rcParams

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['glue-preprocess.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all.h5ad')
push_status(prefix+' read_h5ad')

singlet_cells = pd.read_csv('stats/cells/'+prefix+'-singlets.txt',header=None,sep='\t')

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

adata.obs['is_singlet'] = np.isin(adata.obs['cell'].to_list(),singlet_cells[0].to_list()).tolist()
adata = adata[adata.obs['is_singlet']]
push_status(prefix+' remove_doublets')

if analysis == 'rna':
	adata.layers['counts'] = adata.X.copy()
	sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
	push_status(prefix+' highly_variable_genes')
	
	sc.pp.normalize_total(adata)
	push_status(prefix+' normalize_total')
	
	sc.pp.log1p(adata)
	push_status(prefix+' log1p')
	
	sc.pp.scale(adata)
	push_status(prefix+' scale')
	
	sc.tl.pca(adata, n_comps=100, use_highly_variable=True, svd_solver='auto')
	push_status(prefix+' pca')
else:
	scglue.data.lsi(adata, n_components=100, n_iter=15)
	push_status(prefix+' lsi')

os.makedirs('stats/umap_glue',exist_ok=True)

adata.write_h5ad(filename='hdf5/anndata_glue_'+prefix+'_all_preprocessed.h5ad')
push_status(prefix+' write_h5ad')

# if analysis == 'rna':
# 	sc.pp.neighbors(adata, metric='cosine')
# 	push_status(prefix+' neighbors')
# else:
# 	sc.pp.neighbors(adata, metric='cosine', use_rep='X_lsi')
# 	push_status(prefix+' neighbors')
# 
# sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
# umap02 = adata.obsm['X_umap'].copy()
# np.savetxt('stats/umap_glue/umap02_'+prefix+'_all.txt.gz',umap02,delimiter='\t')
# push_status(prefix+' umap')
# 
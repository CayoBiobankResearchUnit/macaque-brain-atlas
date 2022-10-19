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
# arguments = ['glue-recluster-preprocess.py','rna','rna','12']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3])-1

# Read in main data
adata = ad.read_h5ad('hdf5/subintegration/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
push_status(prefix+' read_h5ad '+prefix+' '+str(this_cluster+1))

if analysis == 'rna':
	# Subset to only good animals
	adata = adata[[x in ['3I4','4I3','6J2'] for x in adata.obs['animal_id']]]
	
	adata.layers['counts'] = adata.X.copy()
	
	try:
		sc.pp.highly_variable_genes(adata, n_top_genes=6000, flavor='seurat_v3')
	except ValueError:
		print('highly_variable_genes failed')
		# var_data = sc.pp.calculate_qc_metrics(adata)[1]
		n_cells = sc.pp.calculate_qc_metrics(adata)[1]['n_cells_by_counts'] > 3
		bdata = adata[:,n_cells.to_list()].copy()
		sc.pp.highly_variable_genes(bdata, n_top_genes=6000, flavor='seurat_v3')
		var_join = bdata.var[['id','highly_variable','highly_variable_rank','means','variances','variances_norm']]
		var = adata.var.merge(var_join,on='id',how='left')
		var.index = var['id'].to_list()
		var['highly_variable'] = [(~np.isnan(x) and x) for x in var['highly_variable']]
		adata.var = var
		# sc.pp.highly_variable_genes(adata, n_top_genes=6000, flavor='seurat_v3')
	
	sc.pp.normalize_total(adata)
	
	sc.pp.log1p(adata)
	
	sc.pp.scale(adata)
	
	sc.tl.pca(adata, n_comps=100, use_highly_variable=True, svd_solver='auto')
	push_status(prefix+' pca '+prefix+' '+str(this_cluster+1))
else:
	scglue.data.lsi(adata, n_components=100, n_iter=15)
	push_status(prefix+' lsi '+prefix+' '+str(this_cluster+1))

adata.write_h5ad(filename='hdf5/subintegration/anndata_'+prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')
push_status(prefix+' write_h5ad '+prefix+' '+str(this_cluster+1))

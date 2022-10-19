#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si, scipy.sparse
import anndata as ad
import scanpy as sc
import numpy as np
import gzip as gz
import os
import shutil
import math

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['atac-concat-cell-motif-matrix.py','atac','atac']

prefix = arguments[1]
analysis = 'atac'

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

all_h5ad = []

for i in sample_list:
	sys.stderr.write('Now appending '+str(i)+'\n')
	this = ad.read_h5ad('stats/motifs/'+i+'-'+prefix+'_cell_by_motif-anndata.h5ad')
	this.X = sp.sparse.csc_matrix(this.X)
	all_h5ad.append(this)

adata = ad.concat(all_h5ad,axis=0,join='inner',merge='same')

obs = pd.read_csv('stats/umap_post/meta_'+prefix+'_all.txt.gz',sep='\t',header=0,index_col=0)

if all(adata.obs.index == obs.index):
	adata.obs = obs

adata.write_h5ad('hdf5/anndata_'+prefix+'_cell_by_motif-all.h5ad')

si.mmwrite('stats/motifs/'+prefix+'_cell_by_motif-all',adata.X,field='integer')

adata.obs['cell'].to_csv('stats/motifs/'+prefix+'_cell_by_motif-all-rows.txt.gz',sep='\t',header=None,index=None)
adata.var['motif'].to_csv('stats/motifs/'+prefix+'_cell_by_motif-all-columns.txt.gz',sep='\t',header=None,index=None)
# adata.obs.to_csv('stats/motifs/'+prefix+'_cell_by_motif-all-metadata.txt.gz',sep='\t',header=True,index=True)

# Gzip the file
with open(os.path.join('stats/motifs',prefix+'_cell_by_motif-all.mtx'),'rb') as mtx_in:
	with gz.open(os.path.join('stats/motifs',prefix+'_cell_by_motif-all.mtx') + '.gz','wb') as mtx_gz:
		shutil.copyfileobj(mtx_in, mtx_gz)
# os.remove(os.path.join(prefix,'matrix.mtx'))

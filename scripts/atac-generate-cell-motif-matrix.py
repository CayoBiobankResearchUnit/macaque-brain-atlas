#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import numpy as np
import gzip as gz
import os
import math

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['atac-generate-cell-motif-matrix.py','atac','atac','1']

prefix = arguments[1]
analysis = 'atac'

# Read in the ID column from the sample metadata data frame
sample_list = pd.read_csv('data/'+prefix+'_metadata.txt',sep='\t').id.tolist()

# Calculate the sample ID (NSM number)
query_sample = sample_list[int(arguments[3])-1]

singlet_cells = pd.read_csv('stats/cells/'+prefix+'-singlets.txt',header=None,sep='\t')

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'-'+query_sample+'.h5ad')

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

adata.obs.index = adata.obs['cell'].to_list()

# Keep only singlet cells
adata.obs['is_singlet'] = np.isin(adata.obs['cell'].to_list(),singlet_cells[0].to_list()).tolist()
adata = adata[adata.obs['is_singlet']]

# Read in peak-motif matrix
peak_tf = si.mmread('sparse-matrices/'+prefix+'-peak_tf.mtx.gz')

# Make it an array
peak_tf_array = peak_tf.toarray()

# Read in the motifs and peaks, respectively, that form rows and columns
peak_tf_rows = pd.read_csv('sparse-matrices/'+prefix+'-peak_tf.rows.txt',index_col=None,header=None,sep='\t')
peak_tf_cols = pd.read_csv('sparse-matrices/'+prefix+'-peak_tf.columns.txt',index_col=None,header=None,sep='\t')

# Index the motifs
peak_tf_rows.index = peak_tf_rows[0].to_list()
peak_tf_rows['motif'] = peak_tf_rows[0].to_list()

# Index the peaks and calculate overlap
peak_tf_cols.index = peak_tf_cols[0].to_list()
peak_tf_cols['is_final'] = np.isin(peak_tf_cols.index.to_list(),adata.var.index.to_list()).tolist()

# Reduce the array and columns (must be in that order) to only final peaks
peak_tf_array = peak_tf_array[:,peak_tf_cols['is_final']]
peak_tf_cols = peak_tf_cols[peak_tf_cols['is_final']]

# Make the peaks categorical (for sorting downstream)
peak_tf_cols[0] = pd.Categorical(
	peak_tf_cols[0],
	categories=adata.var.index.to_list()
)

# Sort both the array and the columns by the peaks (force peaks into same order as peak-by-cell matrix)
peak_tf_final = peak_tf_array[:,np.argsort(peak_tf_cols[0])]
peak_tf_cols = peak_tf_cols.iloc[np.argsort(peak_tf_cols[0])]

# Report a warning if orders do not match
if not all(peak_tf_cols.index == adata.var.index):
	sys.exit('Peaks order does not match')

# Calculate the matrix product to make a cell-by-motif matrix
cell_motif_matrix = adata.X.dot(peak_tf_final.transpose())
push_status(prefix+' dot '+query_sample)

os.makedirs('stats/motifs',exist_ok=True)
np.savetxt('stats/motifs/'+query_sample+'-'+prefix+'_cell_by_motif-matrix.txt.gz',cell_motif_matrix.transpose(),fmt='%s',delimiter='\t')

adata.obs['cell'].to_csv('stats/motifs/'+query_sample+'-'+prefix+'_cell_by_motif-columns.txt.gz',sep='\t',header=None,index=None)
peak_tf_rows['motif'].to_csv('stats/motifs/'+query_sample+'-'+prefix+'_cell_by_motif-rows.txt.gz',sep='\t',header=None,index=None)

# motif
mdata = ad.AnnData(cell_motif_matrix,dtype='int',obs=adata.obs[['cell']],var=peak_tf_rows[['motif']])
mdata.write_h5ad(filename='stats/motifs/'+query_sample+'-'+prefix+'_cell_by_motif-anndata.h5ad')

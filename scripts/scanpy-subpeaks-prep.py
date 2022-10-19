#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import numpy as np
import os

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-subpeaks-prep.py','501']

jobid = int(arguments[1])
this_job = pd.read_csv('data/transposition_file_sizes.txt',header=None,sep='\t',index_col=None).iloc[jobid-1]

query_sample = str(this_job[1])
this_cluster = int(this_job[2])-1

# For recycled code
prefix='atac' # could be "biccn" for integrated things. be careful
analysis='atac'

# sparse-matrices/class${cell_class}-matrix/${this}-class${cell_class}-peak_matrix.mtx.gz

# anndata matrices are transposed relative to what we usually use in monocle3
mm = si.mmread('sparse-matrices/class'+str(this_cluster+1)+'-matrix/'+query_sample+'-class'+str(this_cluster+1)+'-peak_matrix.mtx.gz').T.tocsc()

cells = pd.read_csv('sparse-matrices/class'+str(this_cluster+1)+'-matrix/'+query_sample+'-class'+str(this_cluster+1)+'-peak_matrix.columns.txt',index_col=None,header=None,sep='\t')[0].to_list()
meta_features = pd.read_csv('sparse-matrices/class'+str(this_cluster+1)+'-matrix/'+query_sample+'-class'+str(this_cluster+1)+'-peak_matrix.rows.txt',index_col=None,header=None,sep='\t')
meta_features.columns = ['peak']

coord_split = meta_features['peak'].str.split(r'[_]')
meta_features['chrom'] = coord_split.map(lambda x: x[0])
meta_features['chromStart'] = coord_split.map(lambda x: x[1])
meta_features['chromEnd'] = coord_split.map(lambda x: x[2])

meta_features['chrom'] = pd.Categorical(
	meta_features['chrom'],
	categories=[str(x) for x in range(1,21)] + ['X','Y']
)
meta_features.index = meta_features['peak'].to_list()

# meta_cells = pd.read_csv('stats/glue/subintegration/meta_biccn_class'+str(this_cluster+1)+'.txt.gz',sep='\t',index_col=0,header=0)
meta_cells = pd.read_csv('stats/clusters/subintegration/'+prefix+'-class'+str(this_cluster+1)+'-celltype-predictions.txt',sep='\t',index_col=0,header=0)

good_cells = list(np.where(list(np.array(np.sum(mm,axis=1)).reshape(-1,)))[0])

# Subset to only cells with peaks
mm = mm[good_cells]
meta_cells = meta_cells.iloc[good_cells]

celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

meta_cells = meta_cells.rename(columns={'glue_subtype':'oldpeaks_glue_subtype','glue_subtype_confidence':'oldpeaks_glue_subtype_confidence'})

# Create an AnnData object from the sparse matrix
adata = ad.AnnData(mm,dtype='int',obs=meta_cells,var=meta_features)

# # Write in all the relevant slots
# adata.obs_names = pd.Index(cells)
# adata.var_names = pd.Index(features)
# 
# Filter
# if analysis == 'atac':
# 	adata = adata[(adata.obs.umi_binarized>=1000) & (adata.obs.umi_binarized <= 1e5) & (adata.obs.FRIP>=0.3)]

# Make an hdf5 directory if needed
os.makedirs('hdf5/subpeaks',exist_ok=True)

# adata.X = adata.X.astype(np.int32)

# Write the AnnData object to an hdf5 file
adata.write_h5ad(filename='hdf5/subpeaks/anndata_'+prefix+'-class'+str(this_cluster+1)+'-'+query_sample+'.h5ad')

# Write the AnnData object to an hdf5 file
# adata.write_loom(filename='hdf5/anndata_'+prefix+'-'+query_sample+'.loom')

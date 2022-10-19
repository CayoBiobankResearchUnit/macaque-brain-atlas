#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import itertools
import numpy as np
import gzip as gz
import os
import networkx as nx
import seaborn as sns
import scglue
import math
from matplotlib import rcParams
import sklearn

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['glue-transfer-labels.py','biccn','rna','atac','1','24']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
do_matrix = int(arguments[4])-1
ncpu = int(arguments[5])

# adata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_combined.h5ad')
# push_status(prefix+' read_h5ad')
# 
# adata.obsm['X_glue'] = np.loadtxt('stats/glue/glue_'+prefix+'_all.txt.gz',dtype='float32')
# push_status(prefix+' loadtxt')
# 
# rn_celltypes = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses.txt',header=None,sep='\t')
# rn_cellsubtypes = pd.read_csv('stats/subclusters/'+rn_prefix+'-cellsubclusters.txt',header=None,sep='\t')
# at_celltypes = pd.read_csv('stats/clusters/'+at_prefix+'-final-cellclasses.txt',header=None,sep='\t')
# 
# rn_celltypes_levels = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
# rn_cellsubtypes_levels = pd.read_csv('stats/subclusters/'+rn_prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()
# at_celltypes_levels = pd.read_csv('stats/clusters/'+at_prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
# 
# rn_metadata = pd.read_csv('data/'+rn_prefix+'_metadata.txt',header=0,sep='\t')
# at_metadata = pd.read_csv('data/'+at_prefix+'_metadata.txt',header=0,sep='\t')
# 
# adata.obs['modality'] = ''
# 
# adata.obs['modality'][np.isin(adata.obs['id'],rn_metadata['id'])] = 'RNA'
# adata.obs['modality'][np.isin(adata.obs['id'],at_metadata['id'])] = 'ATAC'
# 
# rnadata = adata[adata.obs['modality'] == 'RNA']
# atadata = adata[adata.obs['modality'] == 'ATAC']
# 
# rn_celltypes.index = rn_celltypes[0].to_list()
# at_celltypes.index = at_celltypes[0].to_list()
# 
# # Sort by cells to match annotated datasets
# rn_celltypes[0] = pd.Categorical(rn_celltypes[0],categories=rnadata.obs['cell'].to_list(), ordered=False)
# rn_cellsubtypes[0] = pd.Categorical(rn_cellsubtypes[0],categories=rnadata.obs['cell'].to_list(), ordered=False)
# at_celltypes[0] = pd.Categorical(at_celltypes[0],categories=atadata.obs['cell'].to_list(), ordered=False)
# 
# rn_celltypes = rn_celltypes.sort_values(by=[0])
# rn_cellsubtypes = rn_cellsubtypes.sort_values(by=[0])
# at_celltypes = at_celltypes.sort_values(by=[0])
# 
# rnadata.obs['rnatype'] = pd.Categorical(rn_celltypes[1],categories=rn_celltypes_levels,ordered=False)
# rnadata.obs['rnasubtype'] = pd.Categorical(rn_cellsubtypes[1],categories=rn_cellsubtypes_levels,ordered=False)
# atadata.obs['atactype'] = pd.Categorical(at_celltypes[1],categories=at_celltypes_levels,ordered=False)
# 
# # rnadata.obs['_index'] = rnadata.obs['cell']
# # atadata.obs['_index'] = atadata.obs['cell']
# 
# rnadata.obs.drop('_index',axis=1,inplace=True)
# atadata.obs.drop('_index',axis=1,inplace=True)
# 
# rnadata.obs.index = rnadata.obs['cell'].to_list()
# atadata.obs.index = atadata.obs['cell'].to_list()
# 
# rn_animal_id = rn_metadata[['id','animal_id']]
# at_animal_id = at_metadata[['id','animal_id']]
# 
# # Add animal IDs
# tmp = rnadata.obs.merge(rn_animal_id,how='inner',on='id')
# tmp.index = tmp['cell'].to_list()
# rnadata.obs = tmp
# tmp = atadata.obs.merge(at_animal_id,how='inner',on='id')
# tmp.index = tmp['cell'].to_list()
# atadata.obs = tmp
# 
# rnadata.write_h5ad('hdf5/anndata_glue_'+prefix+'_'+rn_prefix+'.h5ad')
# atadata.write_h5ad('hdf5/anndata_glue_'+prefix+'_'+at_prefix+'.h5ad')

rnadata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_'+rn_prefix+'.h5ad')
atadata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_'+at_prefix+'.h5ad')
push_status(prefix+' read_h5ad')

animal_intersect = np.intersect1d(rnadata.obs['animal_id'],atadata.obs['animal_id']).tolist()

# Reference RNA (only matching animals)
rnaref = rnadata[np.isin(rnadata.obs['animal_id'],animal_intersect).tolist()]
rnaqry = rnadata[(~np.isin(rnadata.obs['animal_id'],animal_intersect)).tolist()]

ataqry = atadata # sc.pp.subsample(atadata,n_obs=100000,random_state=0,copy=True)

n_sample = 100000

# rnaref = sc.pp.subsample(rnaref,n_obs=n_sample,random_state=0,copy=True)

# RNA evaluation dataset 1: subsample 100k cells from reference
rnaqry_in = sc.pp.subsample(rnaref,n_obs=n_sample,random_state=1,copy=True)

# RNA evaluation dataset 2: subsample 100k cells from withheld samples
rnaqry_ex = sc.pp.subsample(rnaqry,n_obs=n_sample,random_state=1,copy=True)

# ncpu = 52 #len(os.sched_getaffinity(0))

rnaref = ad.read_h5ad('hdf5/glue_reference_'+rn_prefix+'.h5ad')

xrep = rnaref.obsm['X_glue']
xnn = sklearn.neighbors.NearestNeighbors(
	n_neighbors=glue_n_neighbors, n_jobs=ncpu
).fit(xrep)

# 0,1,2,3 should be ataqry, 4,5,6 should be rnaqry_ex, 7,8,9 should be rnaqry_in

if math.ceil(do_matrix/3) <= 1:
	query_prefix = at_prefix
elif math.ceil(do_matrix/3) <= 2:
	query_prefix = rn_prefix + 'ext'
elif math.ceil(do_matrix/3) <= 3:
	query_prefix = rn_prefix + 'int'

query = ad.read_h5ad('hdf5/glue_query_'+query_prefix+'.h5ad')

yrep = query.obsm['X_glue']
ynn = sklearn.neighbors.NearestNeighbors(
	n_neighbors=glue_n_neighbors, n_jobs=ncpu
).fit(yrep)

# If 0

if do_matrix == 0:
	xx = xnn.kneighbors_graph(xrep)
	sp.sparse.save_npz('npz/'+rn_prefix+'-glue-xx.npz',xx)
elif do_matrix % 3 == 1:
	xy = ynn.kneighbors_graph(xrep)
	sp.sparse.save_npz('npz/'+query_prefix+'-glue-xy.npz',xy)
elif do_matrix % 3 == 2:
	yx = xnn.kneighbors_graph(yrep)
	sp.sparse.save_npz('npz/'+query_prefix+'-glue-yx.npz',yx)
elif do_matrix % 3 == 0:
	yy = ynn.kneighbors_graph(yrep)
	sp.sparse.save_npz('npz/'+query_prefix+'-glue-yy.npz',yy)

push_status(prefix+' kneighbors_graph ' +str(do_matrix + 1))

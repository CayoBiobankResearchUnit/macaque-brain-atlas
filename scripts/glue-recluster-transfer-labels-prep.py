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
# arguments = ['glue-transfer-labels.py','biccn','rna','atac','11']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
this_cluster = int(arguments[4])-1

adata = ad.read_h5ad('hdf5/subintegration/anndata_glue_'+prefix+'_combined_class'+str(this_cluster+1)+'.h5ad')

adata.obsm['X_glue'] = np.loadtxt('stats/glue/subintegration/glue_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')

rn_celltypes = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses.txt',header=None,sep='\t')
rn_cellsubtypes = pd.read_csv('stats/subclusters/'+rn_prefix+'-cellsubclusters.txt',header=None,sep='\t')
at_celltypes = pd.read_csv('stats/clusters/'+at_prefix+'-final-cellclasses.txt',header=None,sep='\t')

rn_celltypes_levels = pd.read_csv('stats/clusters/'+rn_prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
rn_cellsubtypes_levels = pd.read_csv('stats/subclusters/'+rn_prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()
at_celltypes_levels = pd.read_csv('stats/clusters/'+at_prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

rn_metadata = pd.read_csv('data/'+rn_prefix+'_metadata.txt',header=0,sep='\t')
at_metadata = pd.read_csv('data/'+at_prefix+'_metadata.txt',header=0,sep='\t')

adata.obs['modality'] = ''

adata.obs['modality'][np.isin(adata.obs['id'],rn_metadata['id'])] = 'RNA'
adata.obs['modality'][np.isin(adata.obs['id'],at_metadata['id'])] = 'ATAC'

rnadata = adata[adata.obs['modality'] == 'RNA']
atadata = adata[adata.obs['modality'] == 'ATAC']

rn_celltypes.index = rn_celltypes[0].to_list()
at_celltypes.index = at_celltypes[0].to_list()

# Sort by cells to match annotated datasets
rn_cellsubtypes[0] = pd.Categorical(rn_cellsubtypes[0],categories=rnadata.obs['cell'].to_list(), ordered=False)

rn_cellsubtypes = rn_cellsubtypes.sort_values(by=[0])

rn_cellsubtypes = rn_cellsubtypes[[str(x)!='nan' for x in rn_cellsubtypes[0]]]

rnadata.obs['rnasubtype'] = pd.Categorical(rn_cellsubtypes[1],categories=rn_cellsubtypes_levels,ordered=False)

# rnadata.obs['_index'] = rnadata.obs['cell']
# atadata.obs['_index'] = atadata.obs['cell']

rnadata.obs.index = rnadata.obs['cell'].to_list()
atadata.obs.index = atadata.obs['cell'].to_list()

rnadata.write_h5ad('hdf5/subintegration/anndata_glue_'+prefix+'_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')
atadata.write_h5ad('hdf5/subintegration/anndata_glue_'+prefix+'_'+at_prefix+'_class'+str(this_cluster+1)+'.h5ad')

# rnadata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_'+rn_prefix+'.h5ad')
# atadata = ad.read_h5ad('hdf5/anndata_glue_'+prefix+'_'+at_prefix+'.h5ad')
# push_status(prefix+' read_h5ad')

# Reference RNA (only matching animals)
rnaref = rnadata

ataqry = atadata # sc.pp.subsample(atadata,n_obs=100000,random_state=0,copy=True)

n_sample = 100000

if rnadata.n_obs > n_sample:
	rnaqry = sc.pp.subsample(rnadata,n_obs=n_sample,random_state=1,copy=True)
else:
	rnaqry = rnadata

rnaref.write_h5ad('hdf5/subintegration/glue_'+prefix+'_reference_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')
ataqry.write_h5ad('hdf5/subintegration/glue_'+prefix+'_query_'+at_prefix+'_class'+str(this_cluster+1)+'.h5ad')
rnaqry.write_h5ad('hdf5/subintegration/glue_'+prefix+'_query_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')

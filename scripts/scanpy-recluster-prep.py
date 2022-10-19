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

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['scanpy-preprocess.py','rna','rna','11']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered.h5ad')
push_status(prefix+' read_h5ad cluster_'+str(this_cluster+1))

cell_types = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes.txt',header=None,sep='\t')
cell_types_levels = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes-levels.txt',header=None,sep='\t')[0].to_list()
cell_classes = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
cell_classes_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

cell_clusters = pd.read_csv('stats/clusters/'+prefix+'-final-clusters.txt',header=None,sep='\t')
cell_partitions = pd.read_csv('stats/clusters/'+prefix+'-final-partitions.txt',header=None,sep='\t')

if all(adata.obs.index == cell_clusters[0].to_list()):
	adata.obs['cluster'] = cell_clusters[1].to_list()

if all(adata.obs.index == cell_partitions[0].to_list()):
	adata.obs['partition'] = cell_partitions[1].to_list()

if all(adata.obs.index == cell_types[0].to_list()):
	adata.obs['cell_type'] = pd.Categorical(
		cell_types[1].to_list(),
		categories=cell_types_levels,
		ordered=False)

if all(adata.obs.index == cell_classes[0].to_list()):
	adata.obs['cell_class'] = pd.Categorical(
		cell_classes[1].to_list(),
		categories=cell_classes_levels,
		ordered=False)

# Split by cluster

os.makedirs('hdf5/recluster',exist_ok=True)

this_cell = cell_classes_levels[this_cluster]

bdata = adata[adata.obs.cell_class == this_cell].copy()

bdata.write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
push_status(prefix+' write_h5ad cluster_'+str(this_cluster+1))

os.makedirs('mm/recluster',exist_ok=True)

si.mmwrite(target='mm/recluster/'+prefix+'_class'+str(this_cluster+1)+'.mtx',a=bdata.X)

bdata.obs['cell'].to_csv('mm/recluster/'+prefix+'_class'+str(this_cluster+1)+'_cols.txt.gz',index=False,header=False)
bdata.var['id'].to_csv('mm/recluster/'+prefix+'_class'+str(this_cluster+1)+'_rows.txt.gz',index=False,header=False)

# for i in range(0,len(cell_classes_levels)):
# 	this_cell = cell_classes_levels[i]
# 	# print(this_cell)
# 	this = adata[adata.obs.cell_class == this_cell].write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(i+1)+'.h5ad')


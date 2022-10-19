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
# arguments = ['glue-recluster-prep.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

cell_class_levels = pd.read_csv('stats/clusters/rna-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all.h5ad')
push_status(prefix+' read_h5ad')

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

all_cells = []
for i in range(0,len(cell_class_levels)):
	try:
		cell_types = pd.read_csv('stats/subintegration/'+prefix+'-class'+str(i+1)+'-cells.txt',header=None,sep='\t')
		all_cells.append(cell_types)
	except pd.errors.EmptyDataError:
		print('No cells detected. Skipping '+cell_class_levels[i])

all_cells = pd.concat(all_cells)
all_cells.index = all_cells[0].to_list()
all_cells.columns = ['cell','cell_class']

meta = adata.obs.merge(all_cells,how='left',on='cell',copy=True)
meta.index = meta['cell'].to_list()
adata.obs = meta

adata = adata[[x in cell_class_levels for x in adata.obs['cell_class']]]

os.makedirs('hdf5/subintegration',exist_ok=True)

for i in range(0,len(cell_class_levels)):
	this_cluster = i
	cell_type = cell_class_levels[i]
	bdata = adata[adata.obs['cell_class'] == cell_type].copy()
	if bool(bdata.n_obs):
		bdata.write_h5ad('hdf5/subintegration/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
	push_status(prefix+' write_h5ad '+str(this_cluster+1))

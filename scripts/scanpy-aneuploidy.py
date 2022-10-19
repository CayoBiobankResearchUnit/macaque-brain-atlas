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
import matplotlib.pyplot as plt
import math
import random

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','rna','2C0','STS']

prefix = arguments[1]
analysis = arguments[2]
animal_id = arguments[3]
region = arguments[4]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
push_status(prefix+' read_h5ad')

# clusters = pd.read_csv('stats/clusters/'+prefix+'-final-clusters.txt',header=None,sep='\t')
# celltypes = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes.txt',header=None,sep='\t')
# celltype_levels = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes-levels.txt',header=None,sep='\t')[0].to_list()
# 
# adata.obs['cluster'] = pd.Categorical(
# 	values=clusters[1],
# 	categories=list(np.unique(clusters[1])),
# )
# adata.obs['celltype'] = pd.Categorical(
# 	values=celltypes[1],
# 	categories=celltype_levels,
# )

# Remove unused chromosomes
adata = adata[:,np.isin(adata.var['chromosome'],[str(i) for i in [*range(1,21,1)]+['X']])]

# Make copy of chromosome and integer-ify it
adata.var['chr'] = adata.var['chromosome']
adata.var['chr'].cat.categories = [*range(1,22,1)]
adata.var['chr'] = adata.var['chr'].astype(int)

# Sort genes in order
adata = adata[:,np.lexsort((adata.var['bp1'].to_numpy(),adata.var['chr'].to_numpy()))]

bdata = adata.copy()

# Subset to just one animal
bdata = bdata[bdata.obs.animal_id == animal_id]

# Subset to just one region
bdata = bdata[bdata.obs.region == region]

# Subset to just three cell types
bdata = bdata[np.isin(bdata.obs['celltype'],['excitatory neurons','radial glial cells','mesenchymal stem cells'])]

# Pick high-UMI cells
bdata = bdata[bdata.obs['n.umi'] > 2000]

# Subset to just 20 cells per cell type
out = []
for i in ['excitatory neurons','radial glial cells','mesenchymal stem cells']:
	random.seed(42)
	cell_ids = random.sample(bdata[bdata.obs.celltype == i].obs.index.to_list(),50)
	out = out + cell_ids

bdata = bdata[out]

# Sort cells in order by cell type
bdata = bdata[np.lexsort((bdata.obs['cluster'].to_numpy(),bdata.obs['celltype'].to_numpy()))]

os.makedirs('stats/aneuploidy',exist_ok=True)

si.mmwrite('stats/aneuploidy/'+prefix+'-expression-matrix-'+animal_id+'-'+region+'.mtx',bdata.X)

bdata.obs.to_csv('stats/aneuploidy/'+prefix+'-meta-'+animal_id+'-'+region+'.txt.gz',sep='\t',header=True,index=True)
bdata.var.to_csv('stats/aneuploidy/'+prefix+'-var-'+animal_id+'-'+region+'.txt.gz',sep='\t',header=True,index=True)

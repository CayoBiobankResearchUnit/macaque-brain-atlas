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
# arguments = ['glue-transfer-labels.py','biccn','rna','atac','2']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
do_matrix = int(arguments[4])-1

ref = ad.read_h5ad('hdf5/glue_reference_'+rn_prefix+'.h5ad')

if do_matrix == 0:
	query_prefix = at_prefix
elif do_matrix == 1:
	query_prefix = rn_prefix+'ext'
elif do_matrix == 2:
	query_prefix = rn_prefix+'int'

query = ad.read_h5ad('hdf5/glue_query_'+query_prefix+'.h5ad')

xx = sp.sparse.load_npz('npz/'+rn_prefix+'-glue-xx.npz')
xy = sp.sparse.load_npz('npz/'+query_prefix+'-glue-xy.npz')
yx = sp.sparse.load_npz('npz/'+query_prefix+'-glue-yx.npz')
yy = sp.sparse.load_npz('npz/'+query_prefix+'-glue-yy.npz')

# glue_n_neighbors = 30

n_neighbors = glue_n_neighbors
field = 'rnatype'
key_added = 'glue_type'
use_rep = 'X_glue'

jaccard = (xx @ yx.T) + (xy @ yy.T)
jaccard.data /= 4 * n_neighbors - jaccard.data
normalized_jaccard = jaccard.multiply(1 / jaccard.sum(axis=0))
onehot = sklearn.preprocessing.OneHotEncoder()
xtab = onehot.fit_transform(ref.obs[[field]])
ytab = normalized_jaccard.T @ xtab
pred = pd.Series(
	onehot.categories_[0][ytab.argmax(axis=1).A1],
	index=query.obs_names, dtype=ref.obs[field].dtype
)
conf = pd.Series(
	ytab.max(axis=1).toarray().ravel(),
	index=query.obs_names
)

query.obs[key_added] = pred
query.obs[key_added + "_confidence"] = conf

query.obs.to_csv('stats/clusters/'+query_prefix+'-celltype-predictions.txt',sep='\t',header=True,index=True)

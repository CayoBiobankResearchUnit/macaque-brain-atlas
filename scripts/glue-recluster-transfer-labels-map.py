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
# arguments = ['glue-transfer-labels.py','biccn','rna','atac','11','1']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
this_cluster = int(arguments[4])-1
do_matrix = int(arguments[5])-1

ref = ad.read_h5ad('hdf5/subintegration/glue_'+prefix+'_reference_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')

if do_matrix == 0:
	query_prefix = at_prefix
elif do_matrix == 1:
	query_prefix = rn_prefix

query = ad.read_h5ad('hdf5/subintegration/glue_'+prefix+'_query_'+query_prefix+'_class'+str(this_cluster+1)+'.h5ad')

xx = sp.sparse.load_npz('npz/subintegration/'+prefix+'-'+rn_prefix+'-glue-class'+str(this_cluster+1)+'-xx.npz')
xy = sp.sparse.load_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-xy.npz')
yx = sp.sparse.load_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-yx.npz')
yy = sp.sparse.load_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-yy.npz')

# glue_n_neighbors = 30

n_neighbors = glue_n_neighbors
field = 'rnasubtype'
key_added = 'glue_subtype'
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

os.makedirs('stats/clusters/subintegration',exist_ok=True)

query.obs.to_csv('stats/clusters/subintegration/'+prefix+'-'+query_prefix+'-class'+str(this_cluster+1)+'-celltype-predictions.txt',sep='\t',header=True,index=True)

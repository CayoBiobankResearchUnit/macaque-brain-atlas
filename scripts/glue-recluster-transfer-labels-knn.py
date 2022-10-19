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
# arguments = ['glue-transfer-labels.py','biccn','rna','atac','1','11','24']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
do_matrix = int(arguments[4])-1
this_cluster = int(arguments[5])-1
ncpu = int(arguments[6])

rnaref = ad.read_h5ad('hdf5/subintegration/glue_'+prefix+'_reference_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')
# push_status(prefix+' read_h5ad')

xrep = rnaref.obsm['X_glue']
xnn = sklearn.neighbors.NearestNeighbors(
	n_neighbors=glue_n_neighbors, n_jobs=ncpu
).fit(xrep)

# 0,1,2,3 should be ataqry, 4,5,6 should be rnaqry_ex, 7,8,9 should be rnaqry_in

if math.ceil(do_matrix/3) <= 1:
	query_prefix = at_prefix
elif math.ceil(do_matrix/3) <= 2:
	query_prefix = rn_prefix

query = ad.read_h5ad('hdf5/subintegration/glue_'+prefix+'_query_'+query_prefix+'_class'+str(this_cluster+1)+'.h5ad')

yrep = query.obsm['X_glue']
ynn = sklearn.neighbors.NearestNeighbors(
	n_neighbors=glue_n_neighbors, n_jobs=ncpu
).fit(yrep)

# If 0

os.makedirs('npz/subintegration',exist_ok=True)

if do_matrix == 0:
	xx = xnn.kneighbors_graph(xrep)
	sp.sparse.save_npz('npz/subintegration/'+prefix+'-'+rn_prefix+'-glue-class'+str(this_cluster+1)+'-xx.npz',xx)
elif do_matrix % 3 == 1:
	xy = ynn.kneighbors_graph(xrep)
	sp.sparse.save_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-xy.npz',xy)
elif do_matrix % 3 == 2:
	yx = xnn.kneighbors_graph(yrep)
	sp.sparse.save_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-yx.npz',yx)
elif do_matrix % 3 == 0:
	yy = ynn.kneighbors_graph(yrep)
	sp.sparse.save_npz('npz/subintegration/'+prefix+'-'+query_prefix+'-glue-class'+str(this_cluster+1)+'-yy.npz',yy)

# push_status(prefix+' kneighbors_graph ' +str(do_matrix + 1))

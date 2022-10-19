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
from matplotlib import rcParams

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['glue-model-regulatory-inference.py','biccn','rna','atac']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]

rnadata_var = pd.read_csv('stats/glue/'+prefix+'_'+rn_prefix+'_regulatory_var.txt.gz',sep='\t',header=0,index_col=0)
atadata_var = pd.read_csv('stats/glue/'+prefix+'_'+at_prefix+'_regulatory_var.txt.gz',sep='\t',header=0,index_col=0)

genes = scglue.genomics.Bed(rnadata_var.assign(name=rnadata_var.index.to_list()).query('highly_variable'))
peaks = scglue.genomics.Bed(atadata_var.assign(name=atadata_var.index.to_list()).query('highly_variable'))

tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)

window_size = 150000

if not os.path.exists('stats/glue/'+prefix+'_dist_regulatory.graphml.gz'):
	dist_graph = scglue.genomics.window_graph(
		promoters, peaks, window_size,
		attr_fn=lambda l, r, d: {
			'dist': abs(d),
			'weight': scglue.genomics.dist_power_decay(abs(d)),
			'type': 'dist'
		}
	)
	dist_graph = nx.DiGraph(dist_graph)
	print(dist_graph.number_of_edges())
	
	nx.write_graphml(dist_graph,'stats/glue/'+prefix+'_dist_regulatory.graphml.gz')
else:
	print('Graph found. Loading from file.')
	dist_graph = nx.read_graphml('stats/glue/'+prefix+'_dist_regulatory.graphml.gz',force_multigraph=True)

feature_embeddings = pd.read_csv('stats/glue/'+prefix+'_feature_embeddings_regulatory.txt.gz',sep='\t',header=0, index_col=0)

reg = scglue.genomics.regulatory_inference(
	feature_embeddings.index,
	feature_embeddings.to_numpy(),
	dist_graph.subgraph([*genes.index, *peaks.index]),
	alternative='greater', random_state=0
)
nx.write_graphml(reg,'stats/glue/'+prefix+'_regulatory_regulatory.graphml.gz')
df = nx.to_pandas_edgelist(reg)

df.to_csv('stats/glue/'+prefix+'_regulatory_regulatory.txt.gz',sep='\t',header=True,index=False)

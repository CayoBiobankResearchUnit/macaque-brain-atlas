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
# arguments = ['glue-recluster-model.py','biccn','rna','atac']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]

rnadata = ad.read_h5ad('hdf5/anndata_glue_graph_'+rn_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+rn_prefix)
atadata = ad.read_h5ad('hdf5/anndata_glue_graph_'+at_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+at_prefix)

rnadata.var['chrom'] = rnadata.var['chromosome']
rnadata.var['chromStart'] = rnadata.var['bp1']
rnadata.var['chromEnd'] = rnadata.var['bp2']
rnadata.var['strand'] = rnadata.var['gene_strand']
genes = scglue.genomics.Bed(rnadata.var.assign(name=rnadata.var_names).query('highly_variable'))

coord_split = atadata.var_names.str.split(r'[_]')
atadata.var['chrom'] = coord_split.map(lambda x: x[0])
atadata.var['chromStart'] = coord_split.map(lambda x: x[1])
atadata.var['chromEnd'] = coord_split.map(lambda x: x[2])

peaks = scglue.genomics.Bed(atadata.var.assign(name=atadata.var_names).query('highly_variable'))

tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)

window_size = 150000

if not os.path.exists('stats/glue/'+prefix+'_dist.graphml.gz'):
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
	
	nx.write_graphml(dist_graph,'stats/glue/'+prefix+'_dist.graphml.gz')
else:
	print('Graph found. Loading from file.')
	dist_graph = nx.read_graphml('stats/glue/'+prefix+'_dist.graphml.gz',force_multigraph=True)

feature_embeddings = pd.read_csv('stats/glue/'+prefix+'_feature_embeddings.txt.gz',sep='\t',header=0, index_col=0)

reg = scglue.genomics.regulatory_inference(
	feature_embeddings.index,
	feature_embeddings.to_numpy(),
	dist_graph.subgraph([*genes.index, *peaks.index]),
	alternative='greater', random_state=0
)
nx.write_graphml(reg,'stats/glue/'+prefix+'_regulatory.graphml.gz')
df = nx.to_pandas_edgelist(reg)

df.to_csv('stats/glue/'+prefix+'_regulatory.txt.gz',sep='\t',header=True,index=False)

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
# arguments = ['glue-recluster-model.py','biccn','rna','atac','11']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
this_cluster = int(arguments[4])-1

rnadata = ad.read_h5ad('hdf5/subintegration/anndata_'+rn_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')
atadata = ad.read_h5ad('hdf5/subintegration/anndata_'+at_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')

rnadata.var['chrom'] = rnadata.var['chromosome']
rnadata.var['chromStart'] = rnadata.var['bp1']
rnadata.var['chromEnd'] = rnadata.var['bp2']
rnadata.var['strand'] = rnadata.var['gene_strand']
genes = scglue.genomics.Bed(rnadata.var.assign(name=rnadata.var_names).query('highly_variable'))

# coord_split = atadata.var_names.str.split(r'[:-]')
if not all([i in atadata.var.columns for i in ['chrom','chromStart','chromEnd']]):
	coord_split = atadata.var_names.str.split(r'[_]')
	atadata.var['chrom'] = coord_split.map(lambda x: x[0])
	atadata.var['chromStart'] = coord_split.map(lambda x: x[1])
	atadata.var['chromEnd'] = coord_split.map(lambda x: x[2])

# Run rna_anchored_prior_graph to calculate highly variable ATAC peaks
graph = scglue.genomics.rna_anchored_prior_graph(rnadata, atadata)

peaks = scglue.genomics.Bed(atadata.var.assign(name=atadata.var_names).query('highly_variable'))

tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)

window_size = 150000

if not os.path.exists('stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_dist.graphml.gz'):
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
	
	nx.write_graphml(dist_graph,'stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_dist.graphml.gz')
else:
	print('Graph found. Loading from file.')
	dist_graph = nx.read_graphml('stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_dist.graphml.gz',force_multigraph=True)

# dist = nx.algorithms.bipartite.biadjacency_matrix(dist_graph, genes.index, peaks.index, weight='dist', dtype=np.float32)
# 
# peak_gene_mapping = scglue.genomics.window_graph(peaks, promoters, 0)
# peak_gene_mapping = nx.DiGraph(peak_gene_mapping)
# peak_gene_mapping = nx.to_pandas_edgelist(
#     peak_gene_mapping, source='Peak1', target='Gene'
# ).loc[:, ['Peak1', 'Gene']]

feature_embeddings = pd.read_csv('stats/glue/subintegration/'+prefix+'_feature_embeddings_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=0, index_col=0)

reg = scglue.genomics.regulatory_inference(
	feature_embeddings.index,
	feature_embeddings.to_numpy(),
	dist_graph.subgraph([*genes.index, *peaks.index]),
	alternative='greater', random_state=0
)
nx.write_graphml(reg,'stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_regulatory.graphml.gz')
df = nx.to_pandas_edgelist(reg)

df.to_csv('stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_regulatory.txt.gz',sep='\t',header=True,index=False)

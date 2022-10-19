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
# arguments = ['glue-model.py','biccn','rna','atac']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]

rnadata = ad.read_h5ad('hdf5/anndata_glue_'+rn_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+rn_prefix)
atadata = ad.read_h5ad('hdf5/anndata_glue_'+at_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' read_h5ad '+at_prefix)

scglue.data.get_gene_annotation(
    rnadata, gtf=gtf_file,
    gtf_by='gene_id'
)
push_status(prefix+' get_gene_annotation')

# coord_split = atadata.var_names.str.split(r'[:-]')
coord_split = atadata.var_names.str.split(r'[_]')
atadata.var['chrom'] = coord_split.map(lambda x: x[0])
atadata.var['chromStart'] = coord_split.map(lambda x: x[1])
atadata.var['chromEnd'] = coord_split.map(lambda x: x[2])

graph = scglue.genomics.rna_anchored_prior_graph(rnadata, atadata)
push_status(prefix+' rna_anchored_prior_graph')

# all(graph.has_node(gene) for gene in rnadata.var_names), \
# all(graph.has_node(peak) for peak in atadata.var_names)

os.makedirs('stats/glue',exist_ok=True)

nx.write_graphml(graph, 'stats/glue/'+prefix+'_prior.graphml.gz')
push_status(prefix+' write_graphml')

rnadata.write_h5ad(filename='hdf5/anndata_glue_graph_'+rn_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' write_h5ad '+rn_prefix)
atadata.write_h5ad(filename='hdf5/anndata_glue_graph_'+at_prefix+'_all_preprocessed.h5ad')
push_status(prefix+' write_h5ad '+at_prefix)

# graph = nx.read_graphml('stats/glue/'+prefix+'_prior.graphml.gz',force_multigraph=True)
# push_status(prefix+' read_graphml')

scglue.models.configure_dataset(
    rnadata, 'NB', use_highly_variable=True,
    use_layer='counts', use_rep='X_pca',
    use_batch='sequencing_run_id'
)
push_status(prefix+' configure_dataset '+rn_prefix)

scglue.models.configure_dataset(
    atadata, 'NB', use_highly_variable=True,
    use_rep='X_lsi',
    use_batch='sequencing_run_id'
)
push_status(prefix+' configure_dataset '+at_prefix)

graph = graph.subgraph(itertools.chain(
    rnadata.var.query('highly_variable').index,
    atadata.var.query('highly_variable').index
))
push_status(prefix+' subgraph')

os.makedirs('glue',exist_ok=True)
glue = scglue.models.fit_SCGLUE(
    {'rna': rnadata, 'atac': atadata}, graph,
    fit_kws={'directory': 'glue'}
)
push_status(prefix+' fit_SCGLUE')

os.makedirs('stats/glue',exist_ok=True)
glue.save('stats/glue/'+prefix+'_glue.dill')
# glue = scglue.models.load_model('stats/glue/'+prefix+'_glue.dill')
push_status(prefix+' save.dill')

dx = scglue.models.integration_consistency(
    glue, {'rna': rnadata, 'atac': atadata}, graph,
    count_layers={'rna': 'counts'}
)
dx.to_csv('stats/glue/'+prefix+'_dx_consistency.txt',sep='\t',header=True,index=True)
push_status(prefix+' glue integration_consistency')

rnadata.obsm['X_glue'] = glue.encode_data('rna', rnadata)
atadata.obsm['X_glue'] = glue.encode_data('atac', atadata)
push_status(prefix+' encode_data')

rnadata.obs['modality'] = 'RNA'
atadata.obs['modality'] = 'ATAC'

# Optimize memory before integrating

rna = rnadata[:,0:0].copy()
atac = atadata[:,0:0].copy()

rna.layers = atac.layers

del rnadata
del atadata

adata = ad.concat([rna, atac])
push_status(prefix+' concat')

adata.obs = adata.obs.drop('animal_integer',1)

os.makedirs('stats/glue',exist_ok=True)
adata.obs.to_csv('stats/glue/meta_'+prefix+'_all.txt.gz',sep='\t',header=True,index=True)

# Save cell embeddings
np.savetxt('stats/glue/glue_'+prefix+'_all.txt.gz',adata.obsm['X_glue'],delimiter='\t')

sc.pp.neighbors(adata, use_rep='X_glue', n_neighbors=15, method='umap', metric='cosine')
push_status(prefix+' neighbors 2D')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/glue/umap02_'+prefix+'_all.txt.gz',umap02,delimiter='\t')
push_status(prefix+' umap 2D')

# Redo neighbors and UMAP to make 10-dimensional UMAP
sc.pp.neighbors(adata, use_rep='X_glue', n_neighbors=30, method='umap', metric='cosine')
push_status(prefix+' neighbors 10D')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/glue/umap10_'+prefix+'_all.txt.gz',umap10,delimiter='\t')
push_status(prefix+' umap 10D')

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/anndata_glue_'+prefix+'_combined.h5ad')
push_status(prefix+' write_h5ad final')

# Write feature embeddings
feature_embeddings = glue.encode_graph(graph)
push_status(prefix+' encode_graph')
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.to_csv('stats/glue/'+prefix+'_feature_embeddings.txt.gz',sep='\t',header=True,index=True)
push_status(prefix+' done')

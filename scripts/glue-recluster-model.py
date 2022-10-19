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
# arguments = ['glue-recluster-model.py','biccn','rna','atac','3']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
this_cluster = int(arguments[4])-1

rnadata = ad.read_h5ad('hdf5/subintegration/anndata_'+rn_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')
atadata = ad.read_h5ad('hdf5/subintegration/anndata_'+at_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')

scglue.data.get_gene_annotation(
    rnadata, gtf=gtf_file,
    gtf_by='gene_id'
)

# coord_split = atadata.var_names.str.split(r'[:-]')
if not all([i in atadata.var.columns for i in ['chrom','chromStart','chromEnd']]):
	coord_split = atadata.var_names.str.split(r'[_]')
	atadata.var['chrom'] = coord_split.map(lambda x: x[0])
	atadata.var['chromStart'] = coord_split.map(lambda x: x[1])
	atadata.var['chromEnd'] = coord_split.map(lambda x: x[2])

graph = scglue.genomics.rna_anchored_prior_graph(rnadata, atadata)

# all(graph.has_node(gene) for gene in rnadata.var_names), \
# all(graph.has_node(peak) for peak in atadata.var_names)

os.makedirs('stats/glue/subintegration',exist_ok=True)

nx.write_graphml(graph,'stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_prior.graphml.gz')
# graph = nx.read_graphml('stats/glue/subintegration/'+prefix+'_class'+str(this_cluster+1)+'_prior.graphml.gz',force_multigraph=True)

scglue.models.configure_dataset(
    rnadata, 'NB', use_highly_variable=True,
    use_layer='counts', use_rep='X_pca',
    use_batch='sequencing_run_id'
)

scglue.models.configure_dataset(
    atadata, 'NB', use_highly_variable=True,
    use_rep='X_lsi',
    use_batch='sequencing_run_id'
)

rnadata.var.to_csv('stats/glue/subintegration/var_'+prefix+'_'+rn_prefix+'_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)
atadata.var.to_csv('stats/glue/subintegration/var_'+prefix+'_'+at_prefix+'_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)

graph = graph.subgraph(itertools.chain(
    rnadata.var.query('highly_variable').index,
    atadata.var.query('highly_variable').index
))

os.makedirs('glue/subintegration'+str(this_cluster+1),exist_ok=True)
os.makedirs('stats/glue/subintegration',exist_ok=True)

if not os.path.exists('stats/glue/subintegration/'+prefix+'_glue_class'+str(this_cluster+1)+'.dill'):
	glue = scglue.models.fit_SCGLUE(
		{'rna': rnadata, 'atac': atadata}, graph,
		fit_kws={'directory': 'glue/subintegration'+str(this_cluster+1)}
	)
	push_status(prefix+' fit_SCGLUE '+str(this_cluster+1))
	
	glue.save('stats/glue/subintegration/'+prefix+'_glue_class'+str(this_cluster+1)+'.dill')
else:
	print('Model found. Loading from file.')
	glue = scglue.models.load_model('stats/glue/subintegration/'+prefix+'_glue_class'+str(this_cluster+1)+'.dill')

if not os.path.exists('stats/glue/subintegration/'+prefix+'_dx_consistency_class'+str(this_cluster+1)+'.txt'):
	dx = scglue.models.integration_consistency(
		glue, {'rna': rnadata, 'atac': atadata}, graph,
		count_layers={'rna': 'counts'}
	)
	dx.to_csv('stats/glue/subintegration/'+prefix+'_dx_consistency_class'+str(this_cluster+1)+'.txt',sep='\t',header=True,index=True)
	push_status(prefix+' glue integration_consistency '+str(this_cluster+1))
else:
	print('Integration consistency found. Skipping.')

rnadata.obsm['X_glue'] = glue.encode_data('rna', rnadata)
atadata.obsm['X_glue'] = glue.encode_data('atac', atadata)

rnadata.obs['modality'] = 'RNA'
atadata.obs['modality'] = 'ATAC'

# Optimize memory before integrating

rna = rnadata[:,0:0].copy()
atac = atadata[:,0:0].copy()

rna.layers = atac.layers

del rnadata
del atadata

adata = ad.concat([rna, atac])

if 'animal_integer' in adata.obs.columns:
	adata.obs = adata.obs.drop('animal_integer',1)

# os.makedirs('stats/glue/subintegration',exist_ok=True)
adata.obs.to_csv('stats/glue/subintegration/meta_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)

# Save cell embeddings
np.savetxt('stats/glue/subintegration/glue_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',adata.obsm['X_glue'],delimiter='\t')

sc.pp.neighbors(adata, use_rep='X_glue', n_neighbors=15, method='umap', metric='cosine')

sc.tl.umap(adata, min_dist=0.25, spread=1.0, n_components=2)
umap02 = adata.obsm['X_umap'].copy()
np.savetxt('stats/glue/subintegration/umap02_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',umap02,delimiter='\t')

# Redo neighbors and UMAP to make 10-dimensional UMAP
sc.pp.neighbors(adata, use_rep='X_glue', n_neighbors=30, method='umap', metric='cosine')

sc.tl.umap(adata, min_dist=0, spread=1.0, n_components=10)
umap10 = adata.obsm['X_umap'].copy()
np.savetxt('stats/glue/subintegration/umap10_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',umap10,delimiter='\t')

adata.obsm['X_umap10'] = umap10.copy()
adata.obsm['X_umap'] = umap02.copy()

adata.write_h5ad(filename='hdf5/subintegration/anndata_glue_'+prefix+'_combined_class'+str(this_cluster+1)+'.h5ad')

# Write feature embeddings
feature_embeddings = glue.encode_graph(graph)

feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.to_csv('stats/glue/subintegration/'+prefix+'_feature_embeddings_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)
push_status(prefix+' done '+str(this_cluster+1))

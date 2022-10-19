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
import pyreadr as pyr

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all.h5ad')
push_status(prefix+' read_h5ad')

# Read in cell metadata
obs = pd.read_csv('stats/umap_post/meta_'+prefix+'_all.txt.gz',index_col=0,sep='\t')

# Read in gene metadata
var = pd.read_csv('stats/var_final/var_'+prefix+'_all.txt.gz',index_col=0,sep='\t')

# Read in UMAP
umap = np.loadtxt('stats/umap_post/umap02_'+prefix+'_all.txt.gz',dtype='float32')

# adata_final = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_postprocessed.h5ad')
# push_status(prefix+' read_h5ad')

if pd.Series(adata.obs['cell']).str.startswith('NSM').any():
	sys.stderr.write('Cell names properly formatted.\n')
else:
	sys.stderr.write('Cell names not properly formatted. Appending sample ID.\n')
	adata.obs['cell'] = np.char.add(np.char.add(adata.obs['id'].to_list(),'-'),adata.obs['cell'].to_list())

# Temporary fix
adata.obs['hemisphere'] = pd.Categorical(adata.obs['hemisphere'], categories=['','L','R'], ordered=False)
adata.obs.iloc()[(adata.obs['region']=='CV').to_list(),(adata.obs.columns=='hemisphere')] = ''
adata.obs.iloc()[(adata.obs['region']=='MB').to_list(),(adata.obs.columns=='hemisphere')] = ''

# Keep only cells in final dataset
cells_in = np.isin(adata.obs['cell'].to_list(),obs.index.to_list()).tolist()
genes_in = np.isin(adata.var.index.to_list(),var.index.to_list()).tolist()

adata = adata[cells_in]
# adata = adata[:,genes_in]
push_status(prefix+' filtered')

# Save raw data
adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_filtered.h5ad')
push_status(prefix+' write_h5ad filtered')

# Count normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
push_status(prefix+' normalize_total')

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p')

# Add in the UMAP
adata.obsm['X_umap'] = umap.copy()

# adata = adata[~(adata.obs['manual_doublet'])]

# Rasterize plots

clusters = pd.read_csv('stats/clusters/'+prefix+'-final-clusters.txt',header=None,sep='\t')
celltypes = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes.txt',header=None,sep='\t')
celltype_levels = pd.read_csv('stats/clusters/'+prefix+'-final-celltypes-levels.txt',header=None,sep='\t')[0].to_list()

adata.obs['cluster'] = pd.Categorical(
	values=clusters[1],
	categories=list(np.unique(clusters[1])),
)
adata.obs['celltype'] = pd.Categorical(
	values=celltypes[1],
	categories=celltype_levels,
)


# # Temporary (write an incomplete version ahead of time so plotting marker genes may proceed)
# adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
# push_status(prefix+' write_h5ad markergenes')

sc.tl.rank_genes_groups(adata, groupby='celltype', method='logreg',max_iter=1000,pts=True)
push_status(prefix+' rank_genes_groups')

sc.tl.rank_genes_groups(adata, groupby='celltype', method='wilcoxon',pts=True,key_added='rank_genes_wilcox')
push_status(prefix+' rank_genes_groups')

sc.tl.rank_genes_groups(adata, groupby='celltype', method='t-test',pts=True,key_added='rank_genes_ttest')
push_status(prefix+' rank_genes_groups')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
push_status(prefix+' write_h5ad markergenes')

cell_types = adata.obs.celltype.cat.categories.to_list()

logreg_results = []
for i in range(0,len(cell_types)):
	this = cell_types[i]
	print(this)
	out = []
	for j in range(0,len(adata.uns['rank_genes_groups']['names'])):
		out.append((this,adata.uns['rank_genes_groups']['names'][j][i],adata.uns['rank_genes_groups']['scores'][j][i]))
	logreg_results.append(pd.DataFrame(out,columns=['cell','gene','logreg_score']))

wilcox_results = []
for i in range(0,len(cell_types)):
	this = cell_types[i]
	print(this)
	out = []
	for j in range(0,len(adata.uns['rank_genes_wilcox']['names'])):
		this_gene = adata.uns['rank_genes_wilcox']['names'][j][i]
		out.append((this,this_gene,adata.uns['rank_genes_wilcox']['logfoldchanges'][j][i],adata.uns['rank_genes_wilcox']['pts'][this][this_gene],adata.uns['rank_genes_wilcox']['pts_rest'][this][this_gene],adata.uns['rank_genes_wilcox']['pvals'][j][i],adata.uns['rank_genes_wilcox']['pvals_adj'][j][i],adata.uns['rank_genes_wilcox']['scores'][j][i]))
	wilcox_results.append(pd.DataFrame(out,columns=['cell','gene','logfoldchanges','pts','pts_rest','wilcox_pvals','wilcox_pvals_adj','wilcox_score']))

ttest_results = []
for i in range(0,len(cell_types)):
	this = cell_types[i]
	print(this)
	out = []
	for j in range(0,len(adata.uns['rank_genes_ttest']['names'])):
		this_gene = adata.uns['rank_genes_ttest']['names'][j][i]
		out.append((this,this_gene,adata.uns['rank_genes_ttest']['logfoldchanges'][j][i],adata.uns['rank_genes_ttest']['pts'][this][this_gene],adata.uns['rank_genes_ttest']['pts_rest'][this][this_gene],adata.uns['rank_genes_ttest']['pvals'][j][i],adata.uns['rank_genes_ttest']['pvals_adj'][j][i],adata.uns['rank_genes_ttest']['scores'][j][i]))
	ttest_results.append(pd.DataFrame(out,columns=['cell','gene','logfoldchanges','pts','pts_rest','ttest_pvals','ttest_pvals_adj','ttest_score']))

logreg = pd.concat(logreg_results)
wilcox = pd.concat(wilcox_results)
ttest = pd.concat(ttest_results)

logreg['gene'] = pd.Categorical(logreg['gene'],categories=adata.var.index.to_list(), ordered=False)
wilcox['gene'] = pd.Categorical(wilcox['gene'],categories=adata.var.index.to_list(), ordered=False)
ttest['gene'] = pd.Categorical(ttest['gene'],categories=adata.var.index.to_list(), ordered=False)
logreg['cell'] = pd.Categorical(logreg['cell'],categories=adata.obs['celltype'].cat.categories.to_list(), ordered=False)
wilcox['cell'] = pd.Categorical(wilcox['cell'],categories=adata.obs['celltype'].cat.categories.to_list(), ordered=False)
ttest['cell'] = pd.Categorical(ttest['cell'],categories=adata.obs['celltype'].cat.categories.to_list(), ordered=False)

logreg = logreg.sort_values(by=['cell','gene']).reset_index(drop=True)
wilcox = wilcox.sort_values(by=['cell','gene']).reset_index(drop=True)
ttest = ttest.sort_values(by=['cell','gene']).reset_index(drop=True)


marker_gene_df = pd.concat([
	wilcox,
	ttest[['ttest_pvals','ttest_pvals_adj','ttest_score']],
	logreg['logreg_score']
],axis=1)

# marker_gene_df['logreg_score'] = logreg['logreg_score']

var = adata.var

marker_gene_df = pd.merge(marker_gene_df,adata.var,left_on='gene',right_on='id',sort=False).sort_values(by=['cell','gene']).reset_index(drop=True)

os.makedirs('rds',exist_ok=True)
pyr.write_rds('rds/'+prefix+'_marker_genes.rds',marker_gene_df)


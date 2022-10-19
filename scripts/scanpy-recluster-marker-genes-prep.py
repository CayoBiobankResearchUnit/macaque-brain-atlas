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
# arguments = ['scanpy-preprocess.py','rna','rna','13']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'.h5ad')
push_status(prefix+' read_h5ad cluster_'+str(this_cluster+1))

# Read in UMAP
umap = np.loadtxt('stats/umap_recluster/umap02_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')
pca = np.loadtxt('stats/umap_recluster/pca_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')

# Read in cell annotations
# celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellcluster2s.txt',header=None,sep='\t')
# celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellcluster1s.txt',header=None,sep='\t')
celltypes = pd.read_csv('stats/subclusters/'+prefix+'-class'+str(this_cluster+1)+'-cellclusters.txt',header=None,sep='\t')

# Save a copy of raw matrix
adata.layers['raw'] = adata.X.copy()

# Normalize data
sc.pp.normalize_total(adata, target_sum=1e4) # ,exclude_highly_expressed=True)
push_status(prefix+' normalize_total cluster_'+str(this_cluster+1))

# Save a copy of normalized matrix
adata.layers['normalized'] = adata.X.copy()

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p cluster_'+str(this_cluster+1))

os.makedirs('npy',exist_ok=True)
x_scaled = sc.pp.scale(adata, zero_center=True,copy=True).X
np.save('npy/scaled_'+prefix+'_class'+str(this_cluster+1)+'_markergenes.npy',x_scaled)

# Add in the UMAP
adata.obsm['X_umap'] = umap.copy()
adata.obsm['X_pca'] = pca.copy()

# adata = adata[~(adata.obs['manual_doublet'])]

adata.obs['cell_subcluster'] = pd.Categorical(
	values=celltypes[1]
)

# adata.obs['cell_subcluster']  = adata.obs['cell_type']
cell_types = adata.obs.cell_subcluster.cat.categories.to_list()

# # Temporary (write an incomplete version ahead of time so plotting marker genes may proceed)
# adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
# push_status(prefix+' write_h5ad markergenes')

if len(cell_types) == 1:
	# Make a dummy dataset of just two empty cells (to force rank_genes_groups to calculate pts)
	bdata = adata[0:2].copy()
	bdata.X = sp.sparse.csr_matrix(bdata.X.shape,dtype=bdata.X.dtype)
	bdata.obs['cell_subcluster'] = 'none'
	bdata.obs.index = ['1','2']
	
	bdata = ad.concat([adata,bdata],axis=0,join='inner',merge='same',uns_merge='same')
	bdata.obs['cell_subcluster'] = pd.Categorical(
		values=bdata.obs['cell_subcluster'],
		categories=[cell_types[0],'none']
	)
	sc.tl.rank_genes_groups(bdata, groupby='cell_subcluster', method='wilcoxon',pts=True,key_added='rank_genes_wilcox')
	
	this = cell_types[0]
	out = []
	for i in range(0,len(adata.var.index)):
		out.append((this,bdata.uns['rank_genes_wilcox']['names'][i][0],bdata.uns['rank_genes_wilcox']['logfoldchanges'][i][0],bdata.uns['rank_genes_wilcox']['pts'][this][0],bdata.uns['rank_genes_wilcox']['pts_rest'][this][0],float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')))
	
	marker_gene_df = pd.DataFrame(out,columns=['cell','gene','logfoldchanges','pts','pts_rest','wilcox_pvals','wilcox_pvals_adj','wilcox_score','ttest_pvals','ttest_pvals_adj','ttest_score','logreg_score'])
	
else:
	logreg_results = []
	if len(cell_types) == 2:
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',max_iter=1000,pts=True,key_added='rank_genes_logreg1')
		push_status(prefix+' rank_genes_groups logreg cluster_'+str(this_cluster+1))
		
		# Change category levels
		adata.obs['cell_subcluster'] = pd.Categorical(
			values=adata.obs['cell_subcluster'],
			categories=cell_types[::-1]
		)
		
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',max_iter=1000,pts=True,key_added='rank_genes_logreg2')
		push_status(prefix+' rank_genes_groups logreg cluster_'+str(this_cluster+1))
		
		# Revert category levels
		adata.obs['cell_subcluster'] = pd.Categorical(
			values=adata.obs['cell_subcluster'],
			categories=cell_types
		)
		
		# Add a 'rank_genes_groups' to keep consistent with below
		adata.uns['rank_genes_groups'] = adata.uns['rank_genes_logreg1']
		
		for i in range(0,len(cell_types)):
			this = cell_types[i]
			print(this)
			out = []
			for j in range(0,len(adata.uns['rank_genes_logreg'+str(2-i)]['names'])):
				out.append((this,adata.uns['rank_genes_logreg'+str(2-i)]['names'][j][0],adata.uns['rank_genes_logreg'+str(2-i)]['scores'][j][0]))
			logreg_results.append(pd.DataFrame(out,columns=['cell','gene','logreg_score']))
	else:
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',max_iter=1000,pts=True)
		push_status(prefix+' rank_genes_groups logreg cluster_'+str(this_cluster+1))
		
		for i in range(0,len(cell_types)):
			this = cell_types[i]
			print(this)
			out = []
			for j in range(0,len(adata.uns['rank_genes_groups']['names'])):
				out.append((this,adata.uns['rank_genes_groups']['names'][j][i],adata.uns['rank_genes_groups']['scores'][j][i]))
			logreg_results.append(pd.DataFrame(out,columns=['cell','gene','logreg_score']))
	
	# Wilcoxon and t-tests are fine
	sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='wilcoxon',pts=True,key_added='rank_genes_wilcox')
	push_status(prefix+' rank_genes_groups wilcox cluster_'+str(this_cluster+1))
	
	sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='t-test',pts=True,key_added='rank_genes_ttest')
	push_status(prefix+' rank_genes_groups ttest cluster_'+str(this_cluster+1))
	
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
	logreg['cell'] = pd.Categorical(logreg['cell'],categories=adata.obs['cell_subcluster'].cat.categories.to_list(), ordered=False)
	wilcox['cell'] = pd.Categorical(wilcox['cell'],categories=adata.obs['cell_subcluster'].cat.categories.to_list(), ordered=False)
	ttest['cell'] = pd.Categorical(ttest['cell'],categories=adata.obs['cell_subcluster'].cat.categories.to_list(), ordered=False)
	
	logreg = logreg.sort_values(by=['cell','gene']).reset_index(drop=True)
	wilcox = wilcox.sort_values(by=['cell','gene']).reset_index(drop=True)
	ttest = ttest.sort_values(by=['cell','gene']).reset_index(drop=True)
	
	marker_gene_df = pd.concat([
		wilcox,
		ttest[['ttest_pvals','ttest_pvals_adj','ttest_score']],
		logreg['logreg_score']
	],axis=1)

	# marker_gene_df['logreg_score'] = logreg['logreg_score']

# Calculate dendrogram
if len(adata.obs['cell_subcluster'].cat.categories) > 1:
	sc.tl.dendrogram(adata,groupby='cell_subcluster',use_rep='X_pca')

var = adata.var

marker_gene_df = pd.merge(marker_gene_df,adata.var,left_on='gene',right_on='id',sort=False).sort_values(by=['cell','gene']).reset_index(drop=True)

adata.write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'_markergenes.h5ad')
push_status(prefix+' write_h5ad markergenes')

os.makedirs('rds',exist_ok=True)
pyr.write_rds('rds/'+prefix+'_class'+str(this_cluster+1)+'_marker_genes.rds',marker_gene_df)


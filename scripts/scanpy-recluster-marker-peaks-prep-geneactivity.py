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
# arguments = ['scanpy-preprocess.py','atac','atac','6']

prefix = arguments[1]
analysis = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_geneactivity.h5ad')
push_status(prefix+' read_h5ad cluster_'+str(this_cluster+1))

# Count normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
push_status(prefix+' normalize_total')

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p')

sc.pp.scale(adata, zero_center=True)


# # Read in UMAP
# umap = np.loadtxt('stats/umap_recluster/umap02_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')
# pca = np.loadtxt('stats/umap_recluster/pca_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')

celltypes = pd.read_csv('stats/clusters/subintegration/'+'subpeak'+'-'+prefix+'sub-class'+str(this_cluster+1)+'-cellsubtype-predictions.txt',header=0,sep='\t',index_col=0)

celltypes = celltypes[['modality','glue_subtype','glue_subtype_confidence']]

celltypes['cell'] = celltypes.index.to_list()

# celltypes = celltypes[celltypes['glue_subtype_confidence'] >= 0.75]

tmp = adata.obs.merge(celltypes,how='left',on='cell',copy=True)
tmp.index = tmp['cell'].to_list()

adata.obs = tmp

cell_subtypes = pd.read_csv('stats/subclusters/'+'rna-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

# adata = adata[adata.obs['glue_subtype_confidence'] >= 0.75]

adata.obs['cell_subcluster'] = pd.Categorical(
	values = adata.obs['glue_subtype'],
	categories = cell_subtypes
)

adata.obs['cell_subcluster'] = adata.obs['cell_subcluster'].cat.remove_unused_categories()

# sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='wilcoxon',pts=True,key_added='rank_peaks_wilcox')

cell_subtypes = adata.obs['cell_subcluster'].cat.categories.to_list()

good_cell_threshold = 50

cell_types = cell_subtypes

adata_backup = adata.copy()
adata = adata[adata.obs['modality'] == 'ATAC']

if any(adata.obs.groupby('cell_subcluster').size() < good_cell_threshold):
	cell_counts = adata.obs.groupby('cell_subcluster').size()
	bad_cells = cell_counts[cell_counts < good_cell_threshold].index.to_list()
	good_cells = cell_counts[cell_counts >= good_cell_threshold].index.to_list()
else:
	good_cells = cell_subtypes
	bad_cells = []

# # Temporary (write an incomplete version ahead of time so plotting marker genes may proceed)
# adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
# push_status(prefix+' write_h5ad markergenes')



adata_2 = adata.copy()

adata_2 = adata_2[pd.Series([x in good_cells for x in adata_2.obs['cell_subcluster']])]
adata_2.obs['purkinje'] = pd.Categorical(adata_2.obs['cell_subcluster'] == 'cerebellar neurons 16')

adata_2.obs['purkinje'] = adata_2.obs['purkinje'].cat.rename_categories({True:'Purkinje',False:'not Purkinje'})

sc.tl.rank_genes_groups(adata_2, groupby='purkinje', method='logreg',groups=['Purkinje','not Purkinje'], reference='rest', max_iter=1000,pts=True,key_added='rank_genes_logreg',multi_class='ovr',solver='lbfgs',n_jobs=-1)



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
	sc.tl.rank_genes_groups(bdata, groupby='cell_subcluster', groups=good_cells, reference='rest'  method='wilcoxon',pts=True,key_added='rank_genes_wilcox')
	
	this = cell_types[0]
	out = []
	for i in range(0,len(adata.var.index)):
		out.append((this,bdata.uns['rank_genes_wilcox']['names'][i][0],bdata.uns['rank_genes_wilcox']['logfoldchanges'][i][0],bdata.uns['rank_genes_wilcox']['pts'][this][0],bdata.uns['rank_genes_wilcox']['pts_rest'][this][0],float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')))
	
	marker_gene_df = pd.DataFrame(out,columns=['cell','gene','logfoldchanges','pts','pts_rest','wilcox_pvals','wilcox_pvals_adj','wilcox_score','ttest_pvals','ttest_pvals_adj','ttest_score','logreg_score'])
	
else:
	logreg_results = []
	if len(cell_types) == 2:
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',groups=good_cells, reference='rest', max_iter=1000,pts=True,key_added='rank_genes_logreg1')
		push_status(prefix+' rank_genes_groups logreg cluster_'+str(this_cluster+1))
		
		# Change category levels
		adata.obs['cell_subcluster'] = pd.Categorical(
			values=adata.obs['cell_subcluster'],
			categories=cell_types[::-1]
		)
		
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',groups=good_cells, reference='rest', max_iter=1000,pts=True,key_added='rank_genes_logreg2')
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
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='logreg',groups=good_cells, reference='rest', max_iter=1000,pts=True,multi_class='ovr',solver='lbfgs',n_jobs=-1)
		push_status(prefix+' rank_genes_groups logreg cluster_'+str(this_cluster+1))
		
		for i in range(0,len(adata.uns['rank_genes_groups']['names'].dtype.names)):
			this = adata.uns['rank_genes_groups']['names'].dtype.names[i]
			print(this)
			out = []
			for j in range(0,len(adata.uns['rank_genes_groups']['names'])):
				out.append((this,adata.uns['rank_genes_groups']['names'][j][i],adata.uns['rank_genes_groups']['scores'][j][i]))
			logreg_results.append(pd.DataFrame(out,columns=['cell','gene','logreg_score']))
	
	# Wilcoxon and t-tests are fine
	sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='wilcoxon',groups=good_cells, reference='rest', pts=True,key_added='rank_genes_wilcox')
	push_status(prefix+' rank_genes_groups wilcox cluster_'+str(this_cluster+1))
	
	sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', method='t-test',groups=good_cells, reference='rest', pts=True,key_added='rank_genes_ttest')
	push_status(prefix+' rank_genes_groups ttest cluster_'+str(this_cluster+1))
	
	ttest = sc.get.rank_genes_groups_df(adata,group=None,key='rank_genes_ttest')
	wilcox = sc.get.rank_genes_groups_df(adata,group=None,key='rank_genes_wilcox')
	
	ttest = ttest[['group','names','logfoldchanges','pct_nz_group','pct_nz_reference','pvals','pvals_adj','scores']]
	ttest = ttest.rename(columns={'group':'cell','names':'gene','pct_nz_group':'pts','pct_nz_reference':'pts_rest','pvals':'ttest_pvals','pvals_adj':'ttest_pvals_adj','scores':'ttest_score'}, errors='raise')
	ttest.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_genes_geneactivity_ttest.pkl')
	
	wilcox = wilcox[['group','names','logfoldchanges','pct_nz_group','pct_nz_reference','pvals','pvals_adj','scores']]
	wilcox = wilcox.rename(columns={'group':'cell','names':'gene','pct_nz_group':'pts','pct_nz_reference':'pts_rest','pvals':'wilcox_pvals','pvals_adj':'wilcox_pvals_adj','scores':'wilcox_score'}, errors='raise')
	wilcox.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_genes_geneactivity_wilcox.pkl')
	
	
	
	#
	logreg = pd.concat(logreg_results)
	
	
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

var = adata.var

marker_gene_df = pd.merge(marker_gene_df,adata.var,left_on='gene',right_on='id',sort=False).sort_values(by=['cell','gene']).reset_index(drop=True)

os.makedirs('rds',exist_ok=True)
pyr.write_rds('rds/'+prefix+'_class'+str(this_cluster+1)+'_marker_genes_geneactivity.rds',marker_gene_df)

marker_gene_df.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_genes_geneactivity.pkl')

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
# arguments = ['scanpy-recluster-marker-peaks-prep.py','atacsub','subpeak','2']

prefix = arguments[1]
# analysis = arguments[2]
glue_prefix = arguments[2]
this_cluster = int(arguments[3]) - 1

# Read in main data
adata = ad.read_h5ad('hdf5/subintegration/anndata_'+prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')
# push_status(prefix+' read_h5ad cluster_'+str(this_cluster+1))

# Read in cell annotations
celltypes = pd.read_csv('stats/clusters/subintegration/'+glue_prefix+'-'+prefix+'-class'+str(this_cluster+1)+'-cellsubtype-predictions.txt',header=0,sep='\t',index_col=0)

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

if any(adata.obs.groupby('cell_subcluster').size() < good_cell_threshold):
	cell_counts = adata.obs.groupby('cell_subcluster').size()
	bad_cells = cell_counts[cell_counts < good_cell_threshold].index.to_list()
	good_cells = cell_counts[cell_counts >= good_cell_threshold].index.to_list()
else:
	good_cells = cell_subtypes
	bad_cells = []

if len(good_cells) == 1:
 	# Make a dummy dataset of just two empty cells (to force rank_genes_groups to calculate pts)
 	this = good_cells[0]
 	out = []
 	for i in range(0,len(adata.var.index)):
 		out.append((this,adata.var.index[i],float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')))

 	marker_peak_df = pd.DataFrame(out,columns=['cell','peak','logfoldchanges','pts','pts_rest','ttest_pvals','ttest_pvals_adj','ttest_score','logreg_score'])
else:
	sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', use_raw=None,
		groups=good_cells, reference='rest', n_genes=None,
		rankby_abs=False, key_added='rank_peaks_ttest', copy=False, method='t-test_overestim_var',
		corr_method='benjamini-hochberg',pts=True)

	ttest = sc.get.rank_genes_groups_df(adata,group=None,key='rank_peaks_ttest')

	ttest = ttest[['group','names','logfoldchanges','pct_nz_group','pct_nz_reference','pvals','pvals_adj','scores']]
	ttest = ttest.rename(columns={'group':'cell','names':'peak','pct_nz_group':'pts','pct_nz_reference':'pts_rest','pvals':'ttest_pvals','pvals_adj':'ttest_pvals_adj','scores':'ttest_score'}, errors='raise')

	ttest['peak'] = pd.Categorical(ttest['peak'],categories=adata.var.index.to_list(), ordered=False)
	ttest['cell'] = pd.Categorical(ttest['cell'],categories=cell_subtypes, ordered=False)

	ttest = ttest.sort_values(by=['cell','peak']).reset_index(drop=True)

	os.makedirs('pkl',exist_ok=True)
	os.makedirs('rds',exist_ok=True)
	ttest.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks_ttest.pkl')
	pyr.write_rds('rds/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks_ttest.rds',ttest)

	logreg_results = []
	if len(good_cells) == 2:
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', use_raw=False,
			groups=good_cells, reference='rest', n_genes=None,
			rankby_abs=False, key_added='rank_peaks_logreg1', copy=False, method='logreg',
			corr_method='benjamini-hochberg',pts=True,max_iter=1000)

		# Change category levels
		adata.obs['cell_subcluster'] = pd.Categorical(
			values=adata.obs['cell_subcluster'],
			categories=cell_subtypes[::-1]
		)

		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', use_raw=False,
			groups=good_cells, reference='rest', n_genes=None,
			rankby_abs=False, key_added='rank_peaks_logreg2', copy=False, method='logreg',
			corr_method='benjamini-hochberg',pts=True,max_iter=1000)

		# Revert category levels
		adata.obs['cell_subcluster'] = pd.Categorical(
			values=adata.obs['cell_subcluster'],
			categories=cell_subtypes
		)

		# Add a 'rank_genes_groups' to keep consistent with below
		adata.uns['rank_peaks_logreg'] = adata.uns['rank_peaks_logreg1']

		for i in range(0,2):
			this = cell_subtypes[i]
			print(this)
			out = []
			for j in range(0,len(adata.uns['rank_peaks_logreg'+str(2-i)]['names'])):
				out.append((this,adata.uns['rank_peaks_logreg'+str(2-i)]['names'][j][0],adata.uns['rank_peaks_logreg'+str(2-i)]['scores'][j][0]))
			logreg_results.append(pd.DataFrame(out,columns=['cell','peak','logreg_score']))
	else:
		# Logistic regression
		sc.tl.rank_genes_groups(adata, groupby='cell_subcluster', use_raw=False,
			groups=good_cells, reference='rest', n_genes=None,
			rankby_abs=False, key_added='rank_peaks_logreg', copy=False, method='logreg',
			corr_method='benjamini-hochberg',pts=True,max_iter=1000,
			multi_class='ovr',solver='lbfgs',n_jobs=-1)

		for i in range(0,len(adata.uns['rank_peaks_logreg']['names'].dtype.names)):
			this = adata.uns['rank_peaks_logreg']['names'].dtype.names[i]
			print(this)
			out = []
			for j in range(0,len(adata.uns['rank_peaks_logreg']['names'])):
				out.append((this,adata.uns['rank_peaks_logreg']['names'][j][i],adata.uns['rank_peaks_logreg']['scores'][j][i]))
			logreg_results.append(pd.DataFrame(out,columns=['cell','peak','logreg_score']))

	logreg = pd.concat(logreg_results)

	logreg['peak'] = pd.Categorical(logreg['peak'],categories=adata.var.index.to_list(), ordered=False)
	logreg['cell'] = pd.Categorical(logreg['cell'],categories=adata.obs['cell_subcluster'].cat.categories.to_list(), ordered=False)

	logreg = logreg.sort_values(by=['cell','peak']).reset_index(drop=True)

	logreg.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks_logreg.pkl')
	pyr.write_rds('rds/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks_logreg.rds',logreg)

	marker_peak_df = ttest.merge(logreg,how='left',on=['cell','peak'],copy=True)

# Calculate dendrogram
if len(adata.obs['cell_subcluster'].cat.categories) > 1:
	sc.tl.dendrogram(adata,groupby='cell_subcluster',use_rep='X_lsi')

adata.write_h5ad(filename='hdf5/recluster/anndata_'+prefix+'_class'+str(this_cluster+1)+'_markerpeaks.h5ad')

os.makedirs('pkl',exist_ok=True)
os.makedirs('rds',exist_ok=True)
marker_peak_df.to_pickle('pkl/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks.pkl')
pyr.write_rds('rds/'+prefix+'_class'+str(this_cluster+1)+'_marker_peaks.rds',marker_peak_df)

os.makedirs('stats/var/subintegration',exist_ok=True)

if not all([x in adata.var.columns for x in ['chrom','chromStart','chromEnd']]):
	coord_split = adata.var_names.str.split(r'[_]')
	adata.var['chrom'] = coord_split.map(lambda x: x[0])
	adata.var['chromStart'] = coord_split.map(lambda x: x[1])
	adata.var['chromEnd'] = coord_split.map(lambda x: x[2])

adata.var[['chrom','chromStart','chromEnd']].to_csv('stats/var/subintegration/var_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',sep='\t',header=True,index=True)

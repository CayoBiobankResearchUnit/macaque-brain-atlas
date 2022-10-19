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
import sklearn as sk
import sklearn.decomposition # TruncatedSVD
import sklearn.preprocessing # normalize
import sklearn.feature_extraction.text # TfidfTransformer

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-recluster-marker-peaks-prep.py','atac','atac']

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

# Read in LSI
lsi = np.loadtxt('stats/preintg_post/pca_'+prefix+'_all.txt.gz',dtype='float32')

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

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_filtered.h5ad')
push_status(prefix+' write_h5ad filtered')

# Read in cell annotations
celltypes = pd.read_csv('stats/clusters/'+prefix+'-celltype-predictions.txt',header=0,sep='\t',index_col=0)
celltype_levels = pd.read_csv('stats/clusters/'+'rna'+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

celltypes = celltypes[['modality','glue_type','glue_type_confidence']]

celltypes = celltypes[celltypes['glue_type_confidence'] >= 0.95]

adata.obs = adata.obs.merge(celltypes,how='left',left_on='cell',right_index=True,copy=True)

adata.obs['glue_type'] = pd.Categorical(
	values = adata.obs['glue_type'],
	categories = celltype_levels
)

# adata = adata[adata.obs['glue_type_confidence'] > 0.95]

adata.obs['glue_type'] = adata.obs['glue_type'].cat.remove_unused_categories()

cell_types = adata.obs['glue_type'].cat.categories.to_list()

# Normalize

## TF-IDF normalization (on highly variable)
#tfidf = sk.feature_extraction.text.TfidfTransformer(norm='l1')
#tfidf.fit(adata.X[:,adata.var['highly_variable']])
#tfidf_matrix = tfidf.transform(adata.X[:,adata.var['highly_variable']])

good_cell_threshold = 50

if any(adata.obs.groupby('glue_type').size() < good_cell_threshold):
	cell_counts = adata.obs.groupby('glue_type').size()
	bad_cells = cell_counts[cell_counts < good_cell_threshold].index.to_list()
	good_cells = cell_counts[cell_counts >= good_cell_threshold].index.to_list()
else:
	good_cells = cell_types
	bad_cells = []

if len(good_cells) == 1:
 	# Make a dummy dataset of just two empty cells (to force rank_genes_groups to calculate pts)
 	this = good_cells[0]
 	out = []
 	for i in range(0,len(adata.var.index)):
 		out.append((this,adata.var.index[i],float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')))
 	
 	marker_peak_df = pd.DataFrame(out,columns=['cell','peak','logfoldchanges','pts','pts_rest','ttest_pvals','ttest_pvals_adj','ttest_score'])
else:
	sc.tl.rank_genes_groups(adata, groupby='glue_type', use_raw=None,
		groups=good_cells, reference='rest', n_genes=None,
		rankby_abs=False, key_added='rank_peaks_ttest', copy=False, method='t-test_overestim_var',
		corr_method='benjamini-hochberg',pts=True)
	
	# ttest_results = []
	# for i in range(0,len(good_cells)):
	# 	this = good_cells[i]
	# 	print(this)
	# 	out = []
	# 	for j in range(0,len(adata.uns['rank_peaks_ttest']['names'])):
	# 		this_peak = adata.uns['rank_peaks_ttest']['names'][j][i]
	# 		out.append((this,this_peak,adata.uns['rank_peaks_ttest']['logfoldchanges'][j][i],adata.uns['rank_peaks_ttest']['pts'][this][this_peak],adata.uns['rank_peaks_ttest']['pts_rest'][this][this_peak],adata.uns['rank_peaks_ttest']['pvals'][j][i],adata.uns['rank_peaks_ttest']['pvals_adj'][j][i],adata.uns['rank_peaks_ttest']['scores'][j][i]))
	# 	ttest_results.append(pd.DataFrame(out,columns=['cell','peak','logfoldchanges','pts','pts_rest','ttest_pvals','ttest_pvals_adj','ttest_score']))
	# 
	# ttest = pd.concat(ttest_results)
	
	ttest = sc.get.rank_genes_groups_df(adata,group=None,key='rank_peaks_ttest')
	
	ttest = ttest[['group','names','logfoldchanges','pct_nz_group','pct_nz_reference','pvals','pvals_adj','scores']]
	ttest = ttest.rename(columns={'group':'cell','names':'peak','pct_nz_group':'pts','pct_nz_reference':'pts_rest','pvals':'ttest_pvals','pvals_adj':'ttest_pvals_adj','scores':'ttest_score'}, errors='raise')
	
	ttest['peak'] = pd.Categorical(ttest['peak'],categories=adata.var.index.to_list(), ordered=False)
	ttest['cell'] = pd.Categorical(ttest['cell'],categories=cell_types, ordered=False)
	
	ttest = ttest.sort_values(by=['cell','peak']).reset_index(drop=True)
	
	os.makedirs('pkl',exist_ok=True)
	os.makedirs('rds',exist_ok=True)
	ttest.to_pickle('pkl/'+prefix+'_marker_peaks_ttest.pkl')
	pyr.write_rds('rds/'+prefix+'_marker_peaks_ttest.rds',ttest)
	
	# Logistic regression
	sc.tl.rank_genes_groups(adata, groupby='glue_type', use_raw=False,
		groups=good_cells, reference='rest', n_genes=None,
		rankby_abs=False, key_added='rank_peaks_logreg', copy=False, method='logreg',
		corr_method='benjamini-hochberg',pts=True,max_iter=1000)
	
	logreg_results = []
	for i in range(0,len(cell_types)):
		this = cell_types[i]
		print(this)
		out = []
		for j in range(0,len(adata.uns['rank_peaks_logreg']['names'])):
			out.append((this,adata.uns['rank_peaks_logreg']['names'][j][i],adata.uns['rank_peaks_logreg']['scores'][j][i]))
		logreg_results.append(pd.DataFrame(out,columns=['cell','peak','logreg_score']))
	
	logreg = pd.concat(logreg_results)
	
	logreg['peak'] = pd.Categorical(logreg['peak'],categories=adata.var.index.to_list(), ordered=False)
	logreg['cell'] = pd.Categorical(logreg['cell'],categories=adata.obs['glue_type'].cat.categories.to_list(), ordered=False)
	
	logreg = logreg.sort_values(by=['cell','peak']).reset_index(drop=True)
	
	logreg.to_pickle('pkl/'+prefix+'_marker_peaks_logreg.pkl')
	pyr.write_rds('rds/'+prefix+'_marker_peaks_logreg.rds',logreg)
	
	marker_peak_df = ttest.merge(logreg,how='left',on=['cell','peak'],copy=True)

# Add in the LSI
adata.obsm['X_lsi'] = lsi.copy()
adata.obsm['X_umap'] = umap.copy()

# Calculate dendrogram
if len(adata.obs['glue_type'].cat.categories) > 1:
	sc.tl.dendrogram(adata,groupby='glue_type',use_rep='X_lsi')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_markerpeaks.h5ad')

os.makedirs('rds',exist_ok=True)
pyr.write_rds('rds/'+prefix+'_marker_peaks.rds',marker_peak_df)

os.makedirs('pkl',exist_ok=True)
marker_peak_df.to_pickle('pkl/'+prefix+'_marker_peaks.pkl')

# Parse peaks to make bed file

coord_split = marker_peak_df['peak'].str.split(r'[_]')
marker_peak_df['chrom'] = pd.Categorical(
	coord_split.map(lambda x: x[0]),
	categories=[str(x) for x in range(1,21)] + ['X','Y']
)
marker_peak_df['chromStart'] = coord_split.map(lambda x: x[1]).astype(int)
marker_peak_df['chromEnd'] = coord_split.map(lambda x: x[2]).astype(int)

# Threshold for calling a significant peak
logreg_threshold = 0.05

os.makedirs('bed',exist_ok=True)
sample_sheet_list = []
for cell in list(marker_peak_df['cell'].unique()):
	print(cell)
	cell_class = celltype_levels.index(cell)+1
	this_cell = marker_peak_df[marker_peak_df['cell'] == cell].sort_values(by=['chrom','chromStart','chromEnd'])
	# this_cell_pos = this_cell[[float(x) > 0 for x in this_cell['ttest_score']]]
	# this_cell_sig = this_cell_pos[this_cell_pos['ttest_pvals_adj'] < enrich_threshold]
	this_cell_sig = this_cell[this_cell['logreg_score'] >= logreg_threshold]
	out_bed = this_cell_sig[['chrom','chromStart','chromEnd','peak','logreg_score']]
	out_file = 'bed/marker_peaks_class'+str(cell_class)+'.bed'
	out_bed.to_csv(out_file,sep='\t',header=None,index=None)
	sample_sheet_list.append([cell,os.getcwd()+'/'+out_file])

sample_sheet = pd.DataFrame(sample_sheet_list).rename(columns={0:'sample_id',1:'sites'})
sample_sheet.to_csv('data/marker_peaks_sample_sheet.tsv',sep='\t',header=True,index=None)

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
# arguments = ['scanpy-marker-genes.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered_endogenous.h5ad')

# Count normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
push_status(prefix+' normalize_total')

# Logarithmize the data
sc.pp.log1p(adata)
push_status(prefix+' log1p')

# adata = adata[~(adata.obs['manual_doublet'])]

# Rasterize plots

cellclasses = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
cellclass_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()

cellclasses[1] = pd.Categorical(
	values=cellclasses[1],
	categories=cellclass_levels,
)

cellclasses = cellclasses.rename(columns={0:'cell',1:'cell_class'})

meta = adata.obs.merge(cellclasses,how='left',on='cell')

missing_cells = meta['cell_class'].value_counts()[meta['cell_class'].value_counts() == 0].index.to_list()

meta['cell_class'] = meta['cell_class'].cat.remove_categories(missing_cells)
meta.index = meta['cell'].to_list()
adata.obs = meta

cellsubclusters = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters.txt',header=None,sep='\t')
cellsubcluster_levels = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

cellsubclusters[1] = pd.Categorical(
	values=cellsubclusters[1],
	categories=cellsubcluster_levels,
)

cellsubclusters = cellsubclusters.rename(columns={0:'cell',1:'cell_subcluster'})

meta = adata.obs.merge(cellsubclusters,how='left',on='cell')

missing_cells = meta['cell_subcluster'].value_counts()[meta['cell_subcluster'].value_counts() == 0].index.to_list()

meta['cell_subcluster'] = meta['cell_subcluster'].cat.remove_categories(missing_cells)
meta.index = meta['cell'].to_list()
adata.obs = meta

push_status(prefix+' prepped')


# # Temporary (write an incomplete version ahead of time so plotting marker genes may proceed)
# adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
# push_status(prefix+' write_h5ad markergenes')

sc.tl.rank_genes_groups(adata, groupby='cell_class', method='logreg',max_iter=1000,pts=True,multi_class='ovr',solver='lbfgs',n_jobs=-1)
push_status(prefix+' rank_genes_groups logreg')

sc.tl.rank_genes_groups(adata, groupby='cell_class', method='wilcoxon',pts=True,key_added='rank_genes_wilcox')
push_status(prefix+' rank_genes_groups wilcox')

sc.tl.rank_genes_groups(adata, groupby='cell_class', method='t-test',pts=True,key_added='rank_genes_ttest')
push_status(prefix+' rank_genes_groups ttest')

adata.write_h5ad(filename='hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
push_status(prefix+' write_h5ad markergenes')

cell_types = adata.obs['cell_class'].cat.categories.to_list()

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
logreg['cell'] = pd.Categorical(logreg['cell'],categories=adata.obs['cell_class'].cat.categories.to_list(), ordered=False)
wilcox['cell'] = pd.Categorical(wilcox['cell'],categories=adata.obs['cell_class'].cat.categories.to_list(), ordered=False)
ttest['cell'] = pd.Categorical(ttest['cell'],categories=adata.obs['cell_class'].cat.categories.to_list(), ordered=False)

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

marker_gene_df.to_pickle('pkl/'+prefix+'_marker_genes.pkl')


# Add in the UMAP
adata.obsm['X_umap'] = np.loadtxt('stats/umap_post/umap02_'+prefix+'_all_endogenous.txt.gz',dtype='float32')

# Read in 10-dimensional UMAP
adata.obsm['X_umap10'] = np.loadtxt('stats/umap_post/umap10_'+prefix+'_all_endogenous.txt.gz',dtype='float32')

# Get a batch-uncorrected PCA
adata.obsm['X_pca'] = np.loadtxt('stats/preintg_post/pca_'+prefix+'_all_endogenous.txt.gz',dtype='float32')


os.makedirs('figures/dendrograms',exist_ok=True)

# sc.pl.dendrogram(adata,groupby='cell_class',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_cellclass.pdf')
# sc.pl.dendrogram(adata,groupby='region',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_region.pdf')
# sc.pl.dendrogram(adata,groupby='cell_subcluster',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_cellsubcluster.pdf')

sce.pp.harmony_integrate(adata,key='sequencing_run_id',basis='X_pca',adjusted_basis='X_pca_harmony')
push_status('harmony integrate')

np.savetxt('stats/umap_post/pca_'+prefix+'_all_harmony_endogenous.txt.gz',adata.obsm['X_pca_harmony'],delimiter='\t')

# adata.obsm['X_pca_harmony'] = np.loadtxt('stats/umap_post/pca_'+prefix+'_all_harmony_endogenous.txt.gz',dtype='float32')


# Function for pulling a scipy linkage object and generating a newick file
def export_newick(anndata, dendrogram_key, groupby, savetree=None):
	# Credit: https://stackoverflow.com/a/31878514
	
	def get_newick(node, parent_dist, leaf_names, newick='') -> str:
		"""
		Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
		
		:param node: output of sciply.cluster.hierarchy.to_tree()
		:param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
		:param leaf_names: list of leaf names
		:param newick: leave empty, this variable is used in recursion.
		:returns: tree in Newick format
		"""
		if node.is_leaf():
			return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
		else:
			if len(newick) > 0:
				newick = "):%.2f%s" % (parent_dist - node.dist, newick)
			else:
				newick = ");"
			newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
			newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
			newick = "(%s" % (newick)
			return newick
	
	tree = scipy.cluster.hierarchy.to_tree(anndata.uns[dendrogram_key]['linkage'],False)
	leaf_names = anndata.obs[groupby].cat.categories.to_list()
	
	if savetree is None:
		return get_newick(tree, tree.dist, leaf_names)
	else:
		with open(savetree, 'w') as f:
			f.write(get_newick(tree, tree.dist, leaf_names)+'\n')


bdata = adata[:,0:0].copy()
# bdata.obsm.clear()
# bdata.obsm['X_pca_harmony'] = adata.obsm['X_pca_harmony']

# Fix region
bdata.obs['region'] = pd.Categorical(
	bdata.obs['region'].astype('str'),
	categories = pd.Index(bdata.obs['region'].cat.categories.to_list() + ['MT']).sort_values().to_list()
)

bdata.obs['region'][(bdata.obs['region'] == 'STS') & (bdata.obs['animal_id'] == '2C0')] = 'MT'

# Fix cell class
keep_cells = bdata.obs['cell_class'].value_counts()[bdata.obs['cell_class'].value_counts() > 2000].index.to_list()
keep_categories = bdata.obs['cell_class'].cat.categories[[i in keep_cells for i in bdata.obs['cell_class'].cat.categories]].to_list()

bdata.obs['cell_class_keep'] = pd.Categorical(
	bdata.obs['cell_class'].astype('str'),
	categories = keep_categories
)

# region_hierarchy = pd.DataFrame({
# 	'region':['dmPFC','vmPFC','dlPFC','vlPFC','ACC','CC','CN','NAc','EC','PC','A1','AMY','HIP','M1','mdTN','vlTN','LGN','S1','IPP','SPP','STS','IT','V1','CV','lCb','MB','MdO','MdC','Pons','MT'],
# 	'region_class':['cortical','cortical','cortical','cortical','cortical','subcortical','subcortical','subcortical','cortical','cortical','cortical','subcortical','subcortical','cortical','subcortical','subcortical','subcortical','cortical','cortical','cortical','cortical','cortical','cortical','cerebellum','cerebellum','brainstem','brainstem','brainstem','brainstem','cortical'],
# 	'region_subclass':['frontal lobe','frontal lobe','frontal lobe','frontal lobe','frontal lobe','corpus callosum','basal ganglia','basal ganglia','temporal lobe','temporal lobe','temporal lobe','hippocampus/amygdala','hippocampus/amygdala','frontal lobe','thalamus','thalamus','thalamus','parietal lobe','parietal lobe','parietal lobe','temporal lobe','temporal lobe','occipital lobe','cerebellum','cerebellum','brainstem','brainstem','brainstem','brainstem','temporal lobe'],
# })
# 
# new_obs = bdata.obs.merge(region_hierarchy,how='outer',on='region')
# new_obs.index = bdata.obs.index
# 
# bdata.obs = new_obs

# bdata.obs['cell_region'] = bdata.obs['cell_class_keep'].astype('str') + '|' + bdata.obs['region_subclass'].astype('str')
# 
# bdata.obs['cell_region'].loc[bdata.obs['cell_region'].str.contains('nan\\|')] = np.NaN
# 
# keep_cell_regions = bdata.obs['cell_region'].value_counts()[bdata.obs['cell_region'].value_counts() > 2000].index.to_list()
# 
# bdata.obs['cell_region_keep'] = pd.Categorical(
# 	bdata.obs['cell_region'].astype('str'),
# 	categories = keep_cell_regions
# )
# 
# bdata.obs['region'] = pd.Categorical(
# 	bdata.obs['region'].astype('str'),
# 	categories = list(np.sort(bdata.obs['region'].unique()))
# )

sc.tl.dendrogram(bdata,groupby='cell_class',use_rep='X_pca_harmony',key_added='dendrogram_pca_cell_class')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_class_keep',use_rep='X_pca_harmony',key_added='dendrogram_pca_cell_class_keep')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_subcluster',use_rep='X_pca_harmony',key_added='dendrogram_pca_cell_subcluster')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='region',use_rep='X_pca_harmony',key_added='dendrogram_pca_region')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_region_keep',use_rep='X_pca_harmony',key_added='dendrogram_pca_cell_region')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_class',use_rep='X_umap10',key_added='dendrogram_umap_cell_class')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_class_keep',use_rep='X_umap10',key_added='dendrogram_umap_cell_class_keep')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_subcluster',use_rep='X_umap10',key_added='dendrogram_umap_cell_subcluster')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='region',use_rep='X_umap10',key_added='dendrogram_umap_region')
push_status(prefix+' dendrogram')

sc.tl.dendrogram(bdata,groupby='cell_region_keep',use_rep='X_umap10',key_added='dendrogram_umap_cell_region')
push_status(prefix+' dendrogram')

# sc.pl.dendrogram(bdata,groupby='cell_class',dendrogram_key='dendrogram_pca_cell_class',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_pca_cellclass.pdf')
# sc.pl.dendrogram(bdata,groupby='cell_class_keep',dendrogram_key='dendrogram_pca_cell_class_keep',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_pca_cellclasskeep.pdf')
# sc.pl.dendrogram(bdata,groupby='region',dendrogram_key='dendrogram_pca_region',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_pca_region.pdf')
# sc.pl.dendrogram(bdata,groupby='cell_subcluster',dendrogram_key='dendrogram_pca_cell_subcluster',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_pca_cellsubcluster.pdf')
# sc.pl.dendrogram(bdata,groupby='cell_region_keep',dendrogram_key='dendrogram_pca_cell_region',orientation='left',save=None).get_figure().savefig('figures/dendrograms/dendrogram_'+prefix+'_pca_cellregion.pdf')


export_newick(bdata,'dendrogram_pca_cell_class_keep','cell_class_keep','stats/dendrograms/'+'dendrogram_pca_cellclasskeep.nwk')
export_newick(bdata,'dendrogram_pca_region','region','stats/dendrograms/'+'dendrogram_pca_region.nwk')

export_newick(bdata,'dendrogram_umap_cell_class_keep','cell_class_keep','stats/dendrograms/'+'dendrogram_umap_cellclasskeep.nwk')
export_newick(bdata,'dendrogram_umap_region','region','stats/dendrograms/'+'dendrogram_umap_region.nwk')

# export_newick(bdata,'dendrogram_pca_cell_region','cell_region_keep')

bdata.write_h5ad(filename='hdf5/anndata_'+prefix+'_slim_dendrogram.h5ad')

# export_newick(bdata,'dendrogram_pca_cell_class','cell_class')
# export_newick(bdata,'dendrogram_pca_region','region')
# 
# for i in range(1,1001):
# 	cdata = bdata[np.random.choice(bdata.n_obs, size=bdata.n_obs, replace=False)]




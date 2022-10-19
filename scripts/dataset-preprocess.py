#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import numpy as np
import gzip as gz
import scrublet as scr
import os
import re

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv

# arguments = ['dataset-preprocess.py','dev-inhibitory-neurons-macaque']
# cortex-dev
# human-cortex
# dev-brain-regions
# adult-brain-vasc-peri
# adult-brain-vasc-endo

prefix = arguments[1]

if prefix in ['cortex-dev']:
	exp_file = 'tsv'
	col_cell = 'Cell'
	col_gene = 'gene'
	col_type = 'WGCNAcluster'
elif prefix in ['human-cortex']:
	exp_file = 'tsv'
	col_cell = 'sample_name'
	col_gene = 'gene'
	col_type = 'cell_type_alias_label'
elif prefix in ['dev-brain-regions']:
	exp_file = 'tsv'
	col_cell = 'cell_name'
	col_gene = 'gene'
	col_type = 'cell_type'
elif prefix in ['adult-brain-vasc-endo']:
	exp_file = 'tsv'
	col_cell = 'Cell'
	col_gene = 'gene'
	col_type = 'Cluster'
elif prefix in ['adult-brain-vasc-peri']:
	exp_file = 'tsv'
	col_cell = 'cellId'
	col_gene = 'gene'
	col_type = 'sub_clusters'
elif prefix in ['vascular-dev']:
	exp_file = 'mtx'
	col_cell = 'cellId'
	col_gene = 'gene'
	col_type = 'cell_type'
elif prefix in ['early-brain']:
	exp_file = 'tsv'
	col_cell = 'Cell'
	col_gene = 'gene'
	col_type = 'Cluster'
elif prefix in ['brain-vasc-atlas']:
	exp_file = 'mtx'
	col_cell = 'Cell'
	col_gene = 'gene'
	col_type = 'Cell_Type'
elif prefix in ['myeloid-neuroinflam']:
	exp_file = 'tsv'
	col_cell = 'ID2'
	col_gene = 'gene'
	col_type = 'Subpopulation'
elif prefix in ['mouse-drg-injury']:
	exp_file = 'tsv'
	col_cell = 'cellId'
	col_gene = 'Feature'
	col_type = 'Classifications'
elif prefix in ['dev-inhibitory-neurons-macaque']:
	exp_file = 'tsv'
	col_cell = 'cellId'
	col_gene = 'gene'
	col_type = 'class'
elif prefix in ['fang-merfish-l2']:
	exp_file = 'tsv'
	col_cell = 'name'
	col_gene = 'name'
	col_type = 'cluster_L2'
elif prefix in ['fang-merfish-l3']:
	exp_file = 'tsv'
	col_cell = 'name'
	col_gene = 'name'
	col_type = 'cluster_L3'

if prefix in ['fang-merfish-l2','fang-merfish-l3']:
	
	adata_list = []
	for i in ['H18.06.006.MTG.4000.expand.rep1','H18.06.006.MTG.4000.expand.rep2','H18.06.006.MTG.4000.expand.rep3','H19.30.001.STG.4000.expand.rep1','H19.30.001.STG.4000.expand.rep2','H20.30.001.STG.4000.expand.rep1','H20.30.001.STG.4000.expand.rep2','H20.30.001.STG.4000.expand.rep3','H22.26.401.MTG.4000.expand.rep1','H22.26.401.MTG.4000.expand.rep2']:
		print(i)
		x = si.mmread('datasets/dryad/'+i+'.matrix.mtx.gz').toarray()
		var = pd.read_csv('datasets/dryad/'+i+'.genes.csv',sep=',',index_col=None,header=0)
		var.index = var[col_gene].to_list()
		meta = pd.read_csv('datasets/dryad/'+i+'.features.csv',index_col=None,sep=',',header=0)
		meta[col_cell] = (i + '-' + meta[col_cell]).to_list()
		meta.index = meta[col_cell].to_list()
		adata = ad.AnnData(x.transpose(),dtype='int',obs=meta,var=var)
		adata_list.append(adata)
	
	adata = ad.concat(adata_list,axis=0,join='inner',merge='same')
	
	# Pseudobulk
	cell_types = list(adata.obs[col_type].unique())
	out_counts = []
	out_cells = []
	for i in cell_types:
		this = adata[adata.obs[col_type] == i]
		out_cells.append(len(this))
		out_counts.append(list(np.sum(this.X,axis=0)))
	
	# Extract gene names
	var = [re.sub('\\|.*','',row) for row in adata.var.index]
	var = pd.DataFrame(var,columns=[col_gene],index=var)
	
	bdata = ad.AnnData(
		np.array(out_counts),
		dtype='int',
		obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells},index=cell_types),
		var = var
	)
	
elif prefix in ['dev-inhibitory-neurons-macaque']:
	adata = ad.read_h5ad('datasets/'+prefix+'.h5ad')
	
	# Pseudobulk
	cell_types = list(adata.obs[col_type].unique())
	out_counts = []
	out_cells = []
	for i in cell_types:
		this = adata[adata.obs[col_type] == i]
		out_cells.append(len(this))
		out_counts.append(list(np.sum(this.raw.X.toarray(),axis=0)))
	
	# Extract gene names
	var = [re.sub('\\|.*','',row) for row in adata.raw.var.index]
	var = pd.DataFrame(var,columns=[col_gene],index=var)
	
	bdata = ad.AnnData(
		np.array(out_counts),
		dtype='int',
		obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells},index=cell_types),
		var = var
	)
elif prefix in ['cortex-dev','human-cortex','dev-brain-regions','adult-brain-vasc-endo','adult-brain-vasc-peri','early-brain','mouse-drg-injury','dev-inhibitory-neurons-macaque','vascular-dev','brain-vasc-atlas','myeloid-neuroinflam']:
	meta = pd.read_csv('datasets/'+prefix+'-meta.tsv',index_col=None,sep='\t')
	meta.index = meta[col_cell].to_list()
	
	# Read in dataset
	if exp_file == 'tsv':
		x = pd.read_csv('datasets/'+prefix+'-expr.tsv.gz',index_col=0,sep='\t')
		
		# Extract gene names
		var = [re.sub('\\|.*','',row) for row in x.index]
		var = pd.DataFrame(var,columns=[col_gene],index=var)
		
		# Convert table to numpy array
		x = x.to_numpy()
	else:
		x = si.mmread('datasets/'+prefix+'-expr.mtx.gz').toarray()
		
		var = pd.read_csv('datasets/'+prefix+'-features.tsv.gz',index_col=None,header=None,sep='\t').rename(columns={0:'gene',1:'gene_name'})
		var.index = var['gene'].to_list()
	
	adata = ad.AnnData(x.transpose(),dtype='int',obs=meta,var=var)
	
	# Pseudobulk
	cell_types = list(adata.obs[col_type].unique())
	out_counts = []
	out_cells = []
	for i in cell_types:
		this = adata[adata.obs[col_type] == i]
		out_cells.append(len(this))
		out_counts.append(list(np.sum(this.X,axis=0)))
	
	bdata = ad.AnnData(
		np.array(out_counts),
		dtype='int',
		obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells},index=cell_types),
		var = var
	)
	
elif prefix == 'rna':
	adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered.h5ad')
	
	# Add cell type annotations
	cell_types = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters.txt',header=None,sep='\t')
	cell_types_levels = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()
	cell_classes = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
	cell_classes_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
	
	cell_classes = cell_classes.rename(columns={0:'cell',1:'cell_class'})
	cell_types = cell_types.rename(columns={0:'cell',1:'cell_subcluster'})
	
	cell_classes['cell_class'] = pd.Categorical(
			cell_classes['cell_class'].to_list(),
			categories=cell_classes_levels,
			ordered=False)
	
	cell_types['cell_subcluster'] = pd.Categorical(
			cell_types['cell_subcluster'].to_list(),
			categories=cell_types_levels,
			ordered=False)
	
	adata.obs = adata.obs.merge(cell_classes,on='cell',how='outer').merge(cell_types,on='cell',how='outer')
	adata.obs.index = adata.obs['cell'].to_list()
	
	adata = adata[~pd.Series(np.isin(adata.obs['cell_class'],['radial glial cells','mesenchymal stem cells','F5 neuron','AHSG neuron','KIR3DL12 neuron','KIR3DL12 microglia']))]
	adata.obs['cell_subcluster'] = adata.obs['cell_subcluster'].cat.remove_unused_categories()
	
	# Pseudobulk (at cell_subcluster level, can also do cell_type, etc.)
	cell_types = adata.obs['cell_subcluster'].cat.categories.to_list()
	out_counts = []
	out_cells = []
	out_classes = []
	out_umi = []
	for i in cell_types:
		this = adata[adata.obs['cell_subcluster'] == i]
		out_classes.append(this.obs['cell_class'].astype('str').unique().tolist()[0])
		out_cells.append(len(this))
		out_counts.append(np.sum(this.X,axis=0).tolist()[0])
		out_umi.append(np.sum(this.X,axis=None))
	
	bdata = ad.AnnData(
		np.array(out_counts),
		dtype='int',
		obs = pd.DataFrame({'celltype':cell_types,'cellclass':out_classes,'n_cells':out_cells,'umi':out_umi},index=cell_types),
		var = adata.var
	)
	
elif prefix == 'rna-class':
	prefix = 'rna'
	
	adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_filtered.h5ad')
	
	# Add cell type annotations
	cell_types = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters.txt',header=None,sep='\t')
	cell_types_levels = pd.read_csv('stats/subclusters/'+prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()
	cell_classes = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses.txt',header=None,sep='\t')
	cell_classes_levels = pd.read_csv('stats/clusters/'+prefix+'-final-cellclasses-levels.txt',header=None,sep='\t')[0].to_list()
	
	cell_classes = cell_classes.rename(columns={0:'cell',1:'cell_class'})
	cell_types = cell_types.rename(columns={0:'cell',1:'cell_subcluster'})
	
	cell_classes['cell_class'] = pd.Categorical(
			cell_classes['cell_class'].to_list(),
			categories=cell_classes_levels,
			ordered=False)
	
	cell_types['cell_subcluster'] = pd.Categorical(
			cell_types['cell_subcluster'].to_list(),
			categories=cell_types_levels,
			ordered=False)
	
	adata.obs = adata.obs.merge(cell_classes,on='cell',how='outer')
	adata.obs.index = adata.obs['cell'].to_list()
	
	adata = adata[~pd.Series(np.isin(adata.obs['cell_class'],['radial glial cells','mesenchymal stem cells','F5 neuron','AHSG neuron','KIR3DL12 neuron','KIR3DL12 microglia']))]
	adata.obs['cell_class'] = adata.obs['cell_class'].cat.remove_unused_categories()
	
	# Pseudobulk (at cell_subcluster level, can also do cell_type, etc.)
	cell_types = adata.obs['cell_class'].cat.categories.to_list()
	out_counts = []
	out_cells = []
	out_classes = []
	out_umi = []
	for i in cell_types:
		this = adata[adata.obs['cell_class'] == i]
		out_cells.append(len(this))
		out_counts.append(np.sum(this.X,axis=0).tolist()[0])
		out_umi.append(np.sum(this.X,axis=None))
	
	bdata = ad.AnnData(
		np.array(out_counts),
		dtype='int',
		obs = pd.DataFrame({'celltype':cell_types,'n_cells':out_cells,'umi':out_umi},index=cell_types),
		var = adata.var
	)
	prefix = 'rna-class'

sc.pp.normalize_total(bdata, target_sum=1e4)
sc.pp.log1p(bdata)

os.makedirs('datasets/processed',exist_ok=True)

np.savetxt('datasets/processed/'+prefix+'-pseudobulk-expr.txt.gz',bdata.X.transpose(),delimiter='\t')
bdata.var.to_csv('datasets/processed/'+prefix+'-pseudobulk-features.txt.gz',sep='\t',header=True,index=True)
bdata.obs.to_csv('datasets/processed/'+prefix+'-pseudobulk-meta.txt.gz',sep='\t',header=True,index=True)

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
import math
from matplotlib import rcParams
import sklearn

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())

arguments = sys.argv
# arguments = ['glue-transfer-labels.py','subpeak','rna','atacsub','15']

prefix = arguments[1]
rn_prefix = arguments[2]
at_prefix = arguments[3]
this_cluster = int(arguments[4])-1

adata = ad.read_h5ad('hdf5/subintegration/anndata_glue_'+prefix+'_combined_class'+str(this_cluster+1)+'.h5ad')

adata.obsm['X_glue'] = np.loadtxt('stats/glue/subintegration/glue_'+prefix+'_class'+str(this_cluster+1)+'.txt.gz',dtype='float32')

rnadata = ad.read_h5ad('hdf5/subintegration/anndata_'+rn_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')
atadata = ad.read_h5ad('hdf5/subintegration/anndata_'+at_prefix+'_class'+str(this_cluster+1)+'_preprocessed.h5ad')

rnaInt = adata[adata.obs['modality'] == 'RNA']
ataInt = adata[adata.obs['modality'] == 'ATAC']

rnadata.obsm['X_glue'] = rnaInt.obsm['X_glue']
atadata.obsm['X_glue'] = ataInt.obsm['X_glue']

if 'modality' in rnadata.obs.columns:
	rnadata.obs = rnadata.obs.drop(['modality'],axis=1)

if 'modality' in atadata.obs.columns:
	atadata.obs = atadata.obs.drop(['modality'],axis=1)

rnaobs = rnadata.obs.merge(rnaInt.obs[['cell']],on='cell',how='left')
rnaobs.index = rnaobs['cell'].to_list()
rnaobs['modality'] = 'RNA'

ataobs = atadata.obs.merge(ataInt.obs[['cell']],on='cell',how='left')
ataobs.index = ataobs['cell'].to_list()
ataobs['modality'] = 'ATAC'

rnadata.obs = rnaobs
atadata.obs = ataobs

# Make counts matrix
rnadata.X = rnadata.layers['counts']

# Add subtype metadata

at_types = pd.read_csv('stats/clusters/subintegration/'+prefix+'-'+at_prefix+'-class'+str(this_cluster+1)+'-cellsubtype-predictions.txt',header=0,sep='\t',index_col=0).rename(columns={'glue_subtype':'cell_subcluster','glue_subtype_confidence':'cell_subcluster_confidence'})
rn_types = pd.read_csv('stats/subclusters/'+rn_prefix+'-class'+str(this_cluster+1)+'-cellclusters.txt',header=None,sep='\t',index_col=0).rename(columns={1:'cell_subcluster'})
rn_types.index = rn_types.index.to_list()

at_types = at_types[['cell_subcluster','cell_subcluster_confidence']]
rn_types['cell_subcluster_confidence'] = np.inf

at_types['cell'] = at_types.index.to_list()
rn_types['cell'] = rn_types.index.to_list()

tmp = rnadata.obs.merge(rn_types,how='left',on='cell',copy=True)
tmp.index = tmp['cell'].to_list()
rnadata.obs = tmp

tmp = atadata.obs.merge(at_types,how='left',on='cell',copy=True)
tmp.index = tmp['cell'].to_list()
atadata.obs = tmp

# celltypes = celltypes[celltypes['glue_subtype_confidence'] >= 0.75]

cell_subtypes = pd.read_csv('stats/subclusters/'+rn_prefix+'-cellsubclusters-levels.txt',header=None,sep='\t')[0].to_list()

# adata = adata[adata.obs['glue_subtype_confidence'] >= 0.75]

rnadata.obs['cell_subcluster'] = pd.Categorical(
	values = rnadata.obs['cell_subcluster'],
	categories = cell_subtypes
).remove_unused_categories()

atadata.obs['cell_subcluster'] = pd.Categorical(
	values = atadata.obs['cell_subcluster'],
	categories = rnadata.obs['cell_subcluster'].cat.categories.to_list()
)

metacell_mean_size = 50
n_meta = math.floor(rnadata.n_obs / metacell_mean_size)

exec(open('scripts/_function_rewrites.py').read())

metacell_results = get_metacells(rnadata,atadata,use_rep='X_glue',n_meta=n_meta,common=True,seed=0)

meta = metacell_results['meta']
adatas = metacell_results['adatas']
meta_summary = []

for i in adatas[0].obs.index.to_list():
	print(i)
	
	this = meta[meta['metacell'] == int(i)]
	
	# pd.crosstab(this['modality'],this['cell_subcluster'])
	
	this_rn = this[this['modality'] == 'RNA']
	this_at = this[this['modality'] == 'ATAC']
	
	rn_counts = this_rn.groupby(['cell_subcluster']).size()
	at_counts = this_at.groupby(['cell_subcluster']).size()
	
	rn_prop = rn_counts / rn_counts.sum()
	at_prop = at_counts / at_counts.sum()
	
	rn_call = this['cell_subcluster'].cat.categories[np.argmax(rn_prop)]
	rn_pct = rn_prop.max()
	rn_score = rn_prop.max() - rn_prop[set(rn_prop.index) ^ set([rn_call])].max()
	
	at_call = this['cell_subcluster'].cat.categories[np.argmax(at_prop)]
	at_pct = at_prop.max()
	at_score = at_prop.max() - at_prop[set(at_prop.index) ^ set([at_call])].max()
	
	m = pd.DataFrame([{
		'metacell_id':str(i),
		'n_rna':this.groupby('modality').size()['RNA'],
		'n_atac':this.groupby('modality').size()['ATAC'],
		'rn_type':rn_call,
		'rn_pct':rn_pct,
		'rn_score':rn_score,
		'at_type':at_call,
		'at_pct':at_pct,
		'at_score':at_score,
		'call_match':rn_call==at_call
	}])
	m.index = m['metacell_id'].to_list()
	
	meta_summary.append(m)

meta_all = pd.concat(meta_summary)

rnameta = adatas[0]
atameta = adatas[1]

rnameta.obs = meta_all
atameta.obs = meta_all

os.makedirs('stats/metacells',exist_ok=True)
os.makedirs('hdf5/metacells',exist_ok=True)

rnameta.write_h5ad('hdf5/metacells/anndata_metacells_'+prefix+'_'+rn_prefix+'_class'+str(this_cluster+1)+'.h5ad')
atameta.write_h5ad('hdf5/metacells/anndata_metacells_'+prefix+'_'+at_prefix+'_class'+str(this_cluster+1)+'.h5ad')

meta_all.to_csv('stats/metacells/metacells_'+prefix+'_class'+str(this_cluster+1)+'_metadata.txt.gz',sep='\t',header=True,index=True)
meta.to_csv('stats/metacells/metacells_'+prefix+'_class'+str(this_cluster+1)+'_metacell_assignments.txt.gz',sep='\t',header=True,index=True)

si.mmwrite('stats/metacells/metacells_'+prefix+'_'+rn_prefix+'_class'+str(this_cluster+1)+'_counts',rnameta.X,field='integer')
si.mmwrite('stats/metacells/metacells_'+prefix+'_'+at_prefix+'_class'+str(this_cluster+1)+'_counts',atameta.X,field='integer')

rnameta.var.to_csv('stats/metacells/metacells_'+prefix+'_'+rn_prefix+'_class'+str(this_cluster+1)+'_var.txt.gz',sep='\t',header=True,index=True)
atameta.var.to_csv('stats/metacells/metacells_'+prefix+'_'+at_prefix+'_class'+str(this_cluster+1)+'_var.txt.gz',sep='\t',header=True,index=True)

# Reduce data

# Temporary for class 1

# if this_cluster == 0:
# 	peaks = pd.read_csv('/scratch/nsnyderm/atac_rna/class1_peaks.txt',index_col=None,header=None,sep='\t')
# 	peaks = pd.DataFrame(list(set(atameta.var.index.to_list()) & set(peaks[0].to_list())))
# 	peaks[0] = pd.Categorical(values=peaks[0],categories=atameta.var.index.to_list()).remove_unused_categories()
# 	peaks = peaks.sort_values(by=0,axis=0)
# 	atameta_sub = atameta[:,peaks[0].to_list()]
# 	si.mmwrite('stats/metacells/metacells_'+at_prefix+'_class'+str(this_cluster+1)+'_counts_sub150k',atameta_sub.X,field='integer')
# 	atameta_sub.var.to_csv('stats/metacells/metacells_'+at_prefix+'_class'+str(this_cluster+1)+'_var_sub150k.txt.gz',sep='\t',header=True,index=True)
	